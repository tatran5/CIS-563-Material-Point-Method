#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class ParticleSystem {
public:
  // TV is a 3D vector of double/float
	using TV = Eigen::Matrix<T,dim,1>;
	using TSM = Eigen::Matrix<T, dim, dim>; // square matrix of values T;

	// Data structures for all particles
  	// Particles have attributes include mass, position, velocity.
	std::vector<T> mParMass;
	std::vector<TV> mParPos;
	std::vector<TV> mParVel;
	T mE;
	T mNu;
	T mMu;
	T mLambda;
	T mRho;
	std::vector<T> mVp0;
	std::vector<TSM> mFp;
	int offsetNode;

  	// Data structures for a Cartesian grid comprising of many nodes. 
  	// (i.e. a cell is marked by 8 nodes (making up a cell (cubde)))
  	// Each node has attribute mass, velocity and position (which is known).
  	// A Node of coordinate (x, y, z) is stored at index
  	// x + y * numCell1D + z * numCell1D * numCell1D
	std::vector<T> mNodeMass;
	std::vector<TV> mNodePos;
	std::vector<TV> mNodeVel;
	
	T mDistBwAdjNode; // distance between two adjacent Nodelines
	int mNumNode1D;
	int mNumParBwAdjNode;


	ParticleSystem() : mE(1e4), mNu(0.3f), mRho(500.f),
			mDistBwAdjNode(0.1f), mNumParBwAdjNode(2), offsetNode(3) {

    	// Setup node attributes
		mNumNode1D = (1.f - 0.f) / mDistBwAdjNode + 1;
		mNodeMass.resize(mNumNode1D * mNumNode1D * mNumNode1D, 0.f);
		mNodePos.resize(mNumNode1D * mNumNode1D * mNumNode1D, TV::Zero());
		mNodeVel.resize(mNumNode1D * mNumNode1D * mNumNode1D, TV::Zero());

	    // Set up particle attributes
	    mMu = mE / (2.f * (1.f + mNu));
	    mLambda = mE * mNu / ((1.f + mNu) * (1.f - 2.f * mNu));

	    // Sample particles randomly between two adjacent Node lines (a cell)
	    // Each cell has 2^dim = 2^3 = 8 particles
    	T approxDistBwPars1D = mDistBwAdjNode / (T) mNumParBwAdjNode; // within each Node
    	// The max offset is the approximate distance between two particles
    	T offsetMax = mDistBwAdjNode / (T) mNumParBwAdjNode;

    	// iterate to the number of Node line - 1 because we do not sample outside of the last Node line
    	for (int xn = offsetNode; xn < mNumNode1D - 1 - offsetNode; xn++) {
    		for (int yn = offsetNode; yn < mNumNode1D - 1 - offsetNode; yn++) {
    			for (int zn = offsetNode; zn < mNumNode1D - 1 - offsetNode; zn++) {
    				// a cell's location is marked by the location of its lower Node lines
    				TV nodeLocWorld = getNodeWorldSpaceFromNodeSpace(xn, yn, zn);			    				

    				for (int xp = 0; xp < mNumParBwAdjNode; xp++) {
    					for (int yp = 0; yp < mNumParBwAdjNode; yp++) {
    						for (int zp = 0; zp < mNumParBwAdjNode; zp++) {
    							mParVel.push_back(TV::Zero());

                				// Setup a particle's location
    							TV parPos = TV::Zero();
    							T offSetX = randTRange(0, offsetMax, 0.000001f);
    							T offSetY = randTRange(0, offsetMax, 0.000001f);
    							T offSetZ = randTRange(0, offsetMax, 0.000001f);

    							parPos(0, 0) = nodeLocWorld(0, 0) + xp * approxDistBwPars1D + offSetX;
    							parPos(1, 0) = nodeLocWorld(1, 0) + yp * approxDistBwPars1D + offSetY;
    							parPos(2, 0) = nodeLocWorld(2, 0) + zp * approxDistBwPars1D + offSetZ;

    							mParPos.push_back(parPos); 							  
    						}
    					}
    				}
    			}
    		}
    	}

    	T numParPerCell = mNumParBwAdjNode * mNumParBwAdjNode * mNumParBwAdjNode;
    	T cellVol = mDistBwAdjNode * mDistBwAdjNode * mDistBwAdjNode;
    	T parVol= cellVol / (T)numParPerCell;

    	int numPar = mParVel.size();
    	mVp0 = std::vector<T>(numPar, parVol);
    	mParMass = std::vector<T>(numPar, parVol * mRho);
    	TSM identityMx = TSM::Identity();
    	mFp = std::vector<TSM>(numPar, identityMx);
    }

    T randTRange(T lo, T hi, T offset) {
    	T loWithOffset = lo + offset;
    	T hiWithOffset = hi - offset;
    	return ((hiWithOffset - loWithOffset) * ((T)std::rand() / RAND_MAX)) + loWithOffset;
    }

    TV getNodeWorldSpaceFromNodeSpace(int x, int y, int z) {
    	TV nodeWorldLoc = TV::Zero();
    	nodeWorldLoc(0, 0) = mDistBwAdjNode * x;
    	nodeWorldLoc(1, 0) = mDistBwAdjNode * y;
    	nodeWorldLoc(2, 0) = mDistBwAdjNode * z;
    	return nodeWorldLoc;
    }

    int getNodeIdxFromNodeSpace(int x, int y, int z) {
    	return x + y * mNumNode1D + z * mNumNode1D * mNumNode1D;
    }


    void computeWeights1D(TV* pWeight, T *pBaseNode, const T& parInNodeSpace) {
    	// Base node is the one closest to the particle in Node space
    	*pBaseNode = std::floor(parInNodeSpace - 0.5) + 1;

    	// Store in the weights that the adjacent nodes in 1D have on this particle
    	T d0, d1, d2, z, z2, zz, zz2;

    	d0 = parInNodeSpace - *pBaseNode + 1;
    	z = 1.5f - d0;
    	z2 = z * z;
    	(*pWeight)(0, 0) = 0.5f * z2;

    	d1 = d0 - 1;
    	(*pWeight)(1, 0) = 0.75f - d1 * d1;

    	d2 = 1.f - d1;
    	zz = 1.5f - d2;
    	zz2 = zz * zz;
    	(*pWeight)(2, 0) = 0.5f * zz2;
    }

    void transferP2G(std::vector<TV>* pActiveNodes, std::vector<TV>* pNodeVelNew) {
    	//std::cout << "P2G 1 " << std::endl;
    	int numPar = mParMass.size();

    	for (int p = 0; p < numPar; p ++) {
    		//std::cout << std::endl;
    		TV parInNodeSpace = mParPos[p] / mDistBwAdjNode;
    		T baseNodeX, baseNodeY, baseNodeZ;
    		TV weightX, weightY, weightZ;
    		computeWeights1D(&weightX, &baseNodeX, parInNodeSpace(0, 0));
    		computeWeights1D(&weightY, &baseNodeY, parInNodeSpace(1, 0));
    		computeWeights1D(&weightZ, &baseNodeZ, parInNodeSpace(2, 0));

    		T sum_w_ijk = 0;
    		for (int i = 0; i < 3; i ++) {
    			T w_i = weightX(i, 0);
    			T node_i = baseNodeX + i;

    			for (int j = 0; j < 3; j++) {
    				T w_ij = w_i * weightY(j, 0);
    				T node_j = baseNodeY + j;
    				
    				for (int k = 0; k < 3; k++) {
    					T w_ijk = w_ij * weightZ(k, 0);
    					T node_k = baseNodeZ + k;
						
						int nodeIdx = getNodeIdxFromNodeSpace(node_i, node_j, node_k);
						
						// splat mass (transfer mass to Node)
	            		mNodeMass[nodeIdx] += mParMass[p] * w_ijk;
    					
	    				// splat momentum (transfer momentum to Node)
						pNodeVelNew->at(nodeIdx) += (w_ijk * mParMass[p]) * mParVel[p];
    					sum_w_ijk += w_ijk;
    				    if (w_ijk < 0) std::cout << "w_ijk" << w_ijk << std::endl; 
                    }
    			}
    		}
    		// Check if sum_w_ijk is 1 -- CORRECT
    		//std::cout << "P2G sum_w_ijk: " << sum_w_ijk << std::endl;
    	}

		//std::cout << "P2G 2 " << std::endl;
    	// update Node velocity by dividing the current momentum by Node mass
    	for (int i = 0; i < mNumNode1D; i++) {
    		for (int j = 0; j < mNumNode1D; j++) {
    			for (int k = 0; k < mNumNode1D; k++) {
    				int nodeIdx = getNodeIdxFromNodeSpace(i, j, k);
    				if (mNodeMass[nodeIdx] != 0) {
    					TV curNodeSpace;
    					curNodeSpace(0, 0) = i;
    					curNodeSpace(1, 0) = j;
    					curNodeSpace(2, 0) = k;
    					pActiveNodes->push_back(curNodeSpace);
    					pNodeVelNew->at(nodeIdx) /= mNodeMass[nodeIdx];
    				} else if (mNodeMass[nodeIdx] < 0) {
                        std::cout << "NEGATIVE NODE MASS IN P2G" << std::endl;
                    } else {
    					pNodeVelNew->at(nodeIdx) = TV::Zero();
    				}
    			}
    		}
    	}
       // std::cout << "P2G 1 " << std::endl;
        
    }

    void computeGravityOnGrid(std::vector<TV>* pForce, const std::vector<TV>& activeNodes, const TV& gravity) {
    //	std::cout << "Grav 1 " << std::endl;
        
        int numActiveNodes = activeNodes.size();
    	for (int i = 0; i < numActiveNodes; i++) {
    		TV activeNode = activeNodes[i];
    		int nodeIdx = getNodeIdxFromNodeSpace(activeNode(0, 0), activeNode(1, 0), activeNode(2, 0));
    		pForce->at(nodeIdx) += mNodeMass[nodeIdx] * gravity;
    	}
      //  std::cout << "Grav 1 " << std::endl;
        
    }

    void polarSVD(TSM* pU, TSM* pSigma, TSM* pV, const TSM& curFp) {
        Eigen::JacobiSVD<TSM> svd(curFp, Eigen::ComputeThinU | Eigen::ComputeThinV);
        *pU = svd.matrixU();
        *pV = svd.matrixV();
        *pSigma = svd.singularValues().asDiagonal();
        
        if (pU->determinant() < 0) {
            (*pU)(0, 2) *= -1;
            (*pU)(1, 2) *= -1;
            (*pU)(2, 2) *= -1;
            (*pSigma)(2, 2) *= -1;
        }

        if (pV->determinant() < 0) {
            (*pV)(0, 2) *= -1;
            (*pV)(1, 2) *= -1;
            (*pV)(2, 2) *= -1;
            (*pSigma)(2, 2) *= -1;
        }
    }

    void fixedCorotated(TSM* pCurP, const TSM& curFp) {
        TSM u, sigma, v;
        polarSVD(&u, &sigma, &v, curFp);
        TSM fChanged = u * sigma * v.transpose();
        TSM r = u * v.transpose();
        T j = fChanged.determinant();
        // std::cout << "u: " << u.determinant() << std::endl;
        // std::cout << "v: " << v.determinant() << std::endl;
        // std::cout << "j: " << j << std::endl;
        TSM a = fChanged.adjoint().transpose();
        *pCurP = 2 * mMu * (fChanged - r) + mLambda * (j - 1) * a;                  
    }

    void computeElasticityOnGrid(std::vector<TV>* pForce) {
       // std::cout << "Ela 1 " << std::endl;
        
        std::vector<TV> before = *pForce;
    	int numPar = mParMass.size();
        for (int p = 0; p < numPar; p++) {
    		TSM curFp = mFp[p];
            TSM curP;
            fixedCorotated(&curP, curFp);
            TSM vp0PFt = mVp0[p] * curP * (curFp).transpose();


            TV curParPos = mParPos[p];
            TV curParPosNodeSpace = curParPos / mDistBwAdjNode;

            TV weightX, weightY, weightZ, dWeightX, dWeightY, dWeightZ;
            T baseNodeX, baseNodeY, baseNodeZ;

            for (int i = 0; i < 3; i++) {
                T wi = weightX(i, 0);
                T dwi_dxi = dWeightX(i, 0) / mDistBwAdjNode;
                T node_i = baseNodeX + i;
                //std::cout << "wi " << wi << std::endl; 

                for (int j = 0; j < 3; j++) {
                    T wj = weightY(j, 0);
                    T wij = wi * wj;
                    T dwij_dxi = dwi_dxi * wj;
                    T dwij_dxj = wi / mDistBwAdjNode * dWeightY(j, 0);
                    T node_j = baseNodeY + j;
                   // std::cout << "wj " << wj << std::endl;

                    for (int k = 0; k < 3; k++) {
                        T wk = weightZ(k, 0);
                        T dwijk_dxi = dwij_dxi * wk;
                        T dwijk_dxj = dwij_dxj * wk;
                        T dwijk_dxk = wij / mDistBwAdjNode * dWeightZ(k, 0);
                        T node_k = baseNodeZ + k;

                        TV gradW = TV(dwijk_dxi, dwijk_dxj, dwijk_dxk);
                        TV foo = -vp0PFt * gradW;
             ///   std::cout << "wk " << wk << std::endl;
                        // std::cout << gradW << std::endl;

                        int nodeIdx = getNodeIdxFromNodeSpace(node_i, node_j, node_k);
                        pForce->at(nodeIdx) += foo;
                    }
                }
            }
    	}
        // for (int i = 0; i < pForce->size(); i ++) {
        //     std::cout << before[i] - pForce->at(i) << std::endl << std::endl;
        // }
      //  std::cout << "Ela 1 " << std::endl;
        
    }

    void setBoundaryVelocities(int thickness) {
        //std::cout << "Bou 1 " << std::endl;
        
    	// send boundary velocity within thickness-range of x direction to be 0 
    	for (int i = 0; i < thickness; i++) {
    		for (int j = 0; j < mNumNode1D; j++) {
    			for (int k = 0; k < mNumNode1D; k++) {
    				int nodeIdx = getNodeIdxFromNodeSpace(i, j, k);
    				mNodeVel[nodeIdx] = TV::Zero();
    			}
    		}
    	} 

    	for (int i = mNumNode1D - thickness; i < mNumNode1D; i++) {
    		for (int j = 0; j < mNumNode1D; j++) {
    			for (int k = 0; k < mNumNode1D; k ++) {
    				int nodeIdx = getNodeIdxFromNodeSpace(i, j, k);
    				mNodeVel[nodeIdx] = TV::Zero();
    			}
    		}
    	} 

    	// send boundary velocity within thickness-range of y direction to be 0 
    	for (int i = 0; i < mNumNode1D; i++) {
    		for (int j = 0; j < thickness; j++) {
    			for (int k = 0; k < mNumNode1D; k++) {
    				int nodeIdx = getNodeIdxFromNodeSpace(i, j, k);
    				mNodeVel[nodeIdx] = TV::Zero();
    			}
    		}
    	} 

    	for (int i = 0; i < mNumNode1D; i ++) {
    		for (int j = mNumNode1D - thickness; j < mNumNode1D; j ++) {
    			for (int k = 0; k < mNumNode1D; k ++) {
    				int nodeIdx = getNodeIdxFromNodeSpace(i, j, k);
    				mNodeVel[nodeIdx] = TV::Zero();
    			}
    		}
    	} 

    	// send boundary velocity within thickness-range of z direction to be 0 
    	for (int i = 0; i < mNumNode1D; i++) {
    		for (int j = 0; j < mNumNode1D; j++) {
    			for (int k = 0; k < thickness; k++) {
    				int nodeIdx = getNodeIdxFromNodeSpace(i, j, k);
    				mNodeVel[nodeIdx] = TV::Zero();
    			}
    		}
    	} 

    	for (int i = 0; i < mNumNode1D; i++) {
    		for (int j = 0; j < mNumNode1D; j++) {
    			for (int k = mNumNode1D - thickness; k < mNumNode1D; k++) {
    				int nodeIdx = getNodeIdxFromNodeSpace(i, j, k);
    				mNodeVel[nodeIdx] = TV::Zero();
    			}
    		}
    	} 
      //  std::cout << "Bou 1 " << std::endl;
        

    }

    void computeWeightsWithGradients1D(TV* pWeight, TV* pDWeight, T *pBaseNode, const T& parInNodeSpace) {
    	// Compute 1D quadratic B spline weights
    	// x is assumed to be scaled in the index space (it is in a grid of size 1)
    	// w is a 3x1 or 1x3 vector

    	*pBaseNode = std::floor(parInNodeSpace - 0.5f);
    	*pWeight= TV::Zero();
    	*pDWeight = TV::Zero();

    	T d0, d1, d2, z, z2, zz, zz2;

    	d0 = parInNodeSpace - *pBaseNode;
    	z = 1.5f - d0;
    	z2 = z * z;
    	(*pWeight)(0, 0) = 0.5f * z2;

    	d1 = d0 - 1.f;
    	(*pWeight)(1, 0) = 0.75f - d1 * d1;

    	d2 = 1 - d1;
    	zz = 1.5f - d2;
    	zz2 = zz * zz;
    	(*pWeight)(2, 0) = 0.5f * zz2;

    	*pDWeight = TV(-z, -2 * d1, zz);
    }

    void evolveF(T dt) {
      //  std::cout << "evolve 1 " << std::endl;
        
    	// evolve deformation gradient
    	int numPar = mParMass.size();
    	TSM curFp, newFp;

    	for (int p = 0; p < numPar; p++) {
    		TV curParPos = mParPos[p];
    		TV curParInNodeSpace = mParPos[p] / mDistBwAdjNode;

    		curFp = mFp[p];

    		T baseNodeX, baseNodeY, baseNodeZ;
    		TV weightX, weightY, weightZ, dWeightX, dWeightY, dWeightZ;

    		computeWeightsWithGradients1D(&weightX, &dWeightX, &baseNodeX, curParInNodeSpace(0, 0));
    		computeWeightsWithGradients1D(&weightY, &dWeightY, &baseNodeY, curParInNodeSpace(1, 0));    		
    		computeWeightsWithGradients1D(&weightZ, &dWeightZ, &baseNodeZ, curParInNodeSpace(2, 0));

    		// Compute grad_vp
    		TSM gradVp = TSM::Zero();
    		for (int i = 0; i < 3; i++) {
    			T wi = weightX(i, 0);
    			T dwi_dxi = dWeightX(i, 0) / mDistBwAdjNode;
    			T node_i = baseNodeX + i;

    			for (int j = 0; j < 3; j++) {
    				T wj = weightY(j, 0);
    				T wij = wi * wj;
    				T dwij_dxi = dwi_dxi * wj;
    				T dwij_dxj = wi / mDistBwAdjNode * dWeightY(j, 0);
    				T node_j = baseNodeY + j;

    				for (int k = 0; k < 3; k++) {
    					T wk = weightZ(k, 0);
    					T wijk = wij * wk;
    					T dwijk_dxi = dwij_dxi * wk;
    					T dwijk_dxj = dwij_dxj * wk;
    					T dwijk_dxk = wij / mDistBwAdjNode * dWeightZ(k, 0);

    					T node_k = baseNodeZ + k;
    					TV gradW = TV(dwijk_dxi, dwijk_dxj, dwijk_dxk);
    					int nodeIdx = getNodeIdxFromNodeSpace(node_i, node_j, node_k);
    					TV v_ijk = mNodeVel[nodeIdx];
    					
    					gradVp += v_ijk * gradW.transpose();
    				}
    			}
    		}

    		newFp = (TSM::Identity() + dt * gradVp) * curFp;

    		for (int i = 0; i < dim; i ++) {
    			for (int j = 0; j < dim; j++) {
    				mFp[p](i, j) = newFp(i, j);
    			}
    		} 
    	}
      //  std::cout << "evolve 1 " << std::endl;
        
    }

    void updateGridVelocity(const std::vector<TV>& activeNodes, const std::vector<TV>& nodeVelNew, 
    	const std::vector<TV>& force, const TV& gravity, const T& dt) {

//std::cout << "update " <<std::endl;

    	int numActiveNodes = activeNodes.size();
    	for (int a = 0; a < numActiveNodes; a++) {
    		int i = activeNodes[a](0, 0);
    		int j = activeNodes[a](1, 0);
    		int k = activeNodes[a](2, 0);
    		int nodeIdx = getNodeIdxFromNodeSpace(i, j, k);
            //std::cout << mNodeMass[nodeIdx] << std::endl;
    		mNodeVel[nodeIdx] = nodeVelNew[nodeIdx] + dt * force[nodeIdx] / mNodeMass[nodeIdx];
    	}
    }

    void transferG2P(const std::vector<TV>& nodeVelNew, const T& flipWeight, const T& dt) {
    	int numPar = mParMass.size();
    	for (int p = 0; p < numPar; p ++) {
    		TV parInNodeSpace = mParPos[p] / mDistBwAdjNode;

    		T baseNodeX, baseNodeY, baseNodeZ;
    		TV weightX, weightY, weightZ;
    		computeWeights1D(&weightX, &baseNodeX, parInNodeSpace(0, 0));
    		computeWeights1D(&weightY, &baseNodeY, parInNodeSpace(1, 0));
    		computeWeights1D(&weightZ, &baseNodeZ, parInNodeSpace(2, 0));	

    		TV vPic = TV::Zero();
    		TV vFlip = mParVel[p];

    		T sum_w_ijk = 0;
    		for (int i = 0; i < 3; i ++) {
    			T w_i = weightX(i, 0);
    			T node_i = baseNodeX + i;

    			for (int j = 0; j < 3; j++) {
    				T w_ij = w_i * weightY(j, 0);
    				T node_j = baseNodeY + j;

    				for (int k = 0; k < 3; k++) {
    					T w_ijk = w_ij * weightZ(k, 0);
    					T node_k = baseNodeZ + k;

    					int nodeIdx = getNodeIdxFromNodeSpace(node_i, node_j, node_k);    					
						vPic += w_ijk * mNodeVel[nodeIdx];
						vFlip += w_ijk * (mNodeVel[nodeIdx] - nodeVelNew[nodeIdx]);

 						sum_w_ijk += w_ijk;
    				}
    			}
    		}

    		// std::cout << "G2P sum_w_ijk: " << (float) sum_w_ijk << std::endl; // debug purpose -- seems to sum up to 1 (desired)

    		mParVel[p] = (1.f - flipWeight) * vPic + flipWeight * vFlip;
    		mParPos[p] += dt * vPic;
    	}
      //  std::cout << "G2P " <<std::endl;
    }

    TV getTotalParticleMomentum() {
    	TV result = TV::Zero();
    	int numPar = mParPos.size();
    	for (int ip = 0; ip < numPar; ip++) {
    		// std::cout << "mParMass[ip]" << mParMass[ip] << std::endl;
    		// std::cout << "mParVel[ip](1, 0)" << mParVel[ip](1, 0) << std::endl;
    		result += mParMass[ip] * mParVel[ip];
    	}
    	return result;
    }

    TV getTotalGridMomentum(const std::vector<T>& nodeMass, const std::vector<TV>& nodeVel) {
    	TV result = TV::Zero();
    	int numNode = nodeMass.size();
    	for (int in = 0; in < numNode; in++) {
    		result += nodeMass[in] * nodeVel[in];
    	}
    	return result;
    }



    void dumpPoly(std::string filename)
    {
    	std::ofstream fs;
    	fs.open(filename);
    	fs << "POINTS\n";
    	int count = 0;
    	for (auto parPos : mParPos) {
    		fs << ++count << ":";
    		for (int i = 0; i < dim; i++){
                if (parPos(i) != parPos(i))
                    std::cout << "GETTING NAN FOR PARTICLE POSITION" << std::endl;
    			fs << " " << parPos(i);
            }
    		if (dim == 2)
    			fs << " 0";
    		fs << "\n";
    	}
    	fs << "POLYS\n";
    	count = 0;
    	fs << "END\n";
    	fs.close();
    }
};
