#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc -- guess thats why not working
#include "tiny_obj_loader.h"
#include "mesh_query0.1/mesh_query.h"

#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include <tbb/tbb.h>


using namespace tinyobj;

template<class T, int dim>
class SimulationDriver{
public:
	using TV = Eigen::Matrix<T,dim,1>;
	using SpMat = Eigen::SparseMatrix<T>;
	using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;
	using TSM = Eigen::Matrix<T, dim, dim>; // square matrix

	// Simulation parameters
	T dt;
	TV gravity;

	// Grid parameters
	T dx; // width of a cell
	T minCorner;
	T maxCorner;
	T res; // number of cells 1D
	// Add-on grid parameters
	T NpPerCell1D;
	T offsetGrid; // so that particles are only sampled within a smaller grid
	T numNode;

	// Particle parameters
	T E; //elasticity
	T nu, mu, lambda, rho, Np, Vp0;
	std::vector<TV> xp; // position of particle
	std::vector<T> mp; // mass of particle
	std::vector<TV> vp; // velocity of particle
	std::vector<TSM> Fp;

	// Other

	SimulationDriver() {
		std::cout << "SimulationDriver called" << std::endl;

		// Set values for simulation parameters
		dt = 1e-3;
		gravity = TV(0, -9.8, 0);

		// Set values for grid parameters
		// dx = 0.02f;
		dx = 0.02f;
		minCorner = 0.f;
		maxCorner = 1.f;
		res = (maxCorner - minCorner) / dx + 1;
		// Set values for add-on grid parameters
		NpPerCell1D = 12.f;
		offsetGrid = 6.f;
		numNode = res * res * res;

		// Set values for particles
		E = 100000;
		nu = 0.3f;
		mu = E / (2 * (1 + nu));
		lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
		rho = 1000;
		//sampleParticles();
		sampleParticlesInMesh();
		Np = xp.size();
		Vp0 = dx * dx * dx / NpPerCell1D / NpPerCell1D / NpPerCell1D;
		mp = std::vector<T>(Np, Vp0 * rho);
		vp = std::vector<TV>(Np, TV::Zero());
		Fp = std::vector<TSM>(Np, TSM::Identity());
		std::cout << "Np = " << Np << std::endl;
		// std::cout << "res = " << res << std::endl;
		// std::cout << "numNode = " << numNode << std::endl;
	}

	void sampleParticles() {
		T estDistBwPars1D = dx / NpPerCell1D;
		for (int cx = offsetGrid; cx < res - 1 - offsetGrid; cx++) {
			for (int cy = offsetGrid; cy < res - 1 - offsetGrid; cy++) {
				for (int cz = offsetGrid; cz < res - 1 - offsetGrid; cz++) {
					TV cellWorldSpace = gridSpace2WorldSpace(cx, cy, cz);
					for (int px = 0; px < NpPerCell1D; px++) {
						for (int py = 0; py < NpPerCell1D; py++) {
							for (int pz = 0; pz < NpPerCell1D; pz++) {
								T offsetX = randWithOffset(0, estDistBwPars1D, 0.000001f);
								T offsetY = randWithOffset(0, estDistBwPars1D, 0.000001f);
								T offsetZ = randWithOffset(0, estDistBwPars1D, 0.000001f);

								TV parPos = TV::Zero();
								parPos[0] = cellWorldSpace[0] + px * estDistBwPars1D + offsetX;
								parPos[1] = cellWorldSpace[1] + py * estDistBwPars1D + offsetY;
								parPos[2] = cellWorldSpace[2] + pz * estDistBwPars1D + offsetZ;
								xp.push_back(parPos);
							}
						}
					}
				}
			}
		}
	}

	void sampleParticlesInMesh() {
		std::cout << "sampleParticlesInMesh called" << std::endl;

		attrib_t attrib;
		std::vector<shape_t> shapes = std::vector<shape_t>();
		std::vector<material_t> materials = std::vector<material_t>();
		std::string warn;
		std::string err;
		std::string filenameStr = "../../Obj/TheFUCK.obj";
		const char* filenameChar = filenameStr.c_str();
		bool isLoaded = LoadObj(&attrib, &shapes, &materials, &warn, &err, filenameChar, NULL, true, true);

		if (!isLoaded) {
			std::cout << "OBJ was not loaded, encountered error: " << err << std::endl;
		} else if (shapes.size() > 1) {
			std::cout << "MORE THAN ONE MESH DETECTED. WE DO NOT ALLOW AN OBJ TO HAVE MORE THAN ONE MESH" << std::endl;
		} else {
			std::cout << "OBJ loaded" << std::endl;
			std::cout << "WARNINGS FROM OBJ LOADER: " << std::endl << warn << "WARNINGS ENDED"<< std::endl;
			// vertices are flatten, so indices 0, 1, 2 holds position x, y, z of the first vertex
			std::vector<real_t>& vertices = attrib.vertices;

			// find the smallest and largest point in each direction
			// so that we can scale this down and put it at position from 0 to 1
			T smallestX = INFINITY;
			T smallestY = INFINITY;
			T smallestZ = INFINITY;
			T largestX = -INFINITY;
			T largestY = -INFINITY;
			T largestZ = -INFINITY;
			for (int v = 0; v < vertices.size(); v++) {
				if ( v % 3 == 0) {
					if (vertices[v] < smallestX) smallestX = vertices[v];
					if (vertices[v] > largestX) largestX = vertices[v];
				} else if ( v % 3 == 1) {
					if (vertices[v] < smallestY) smallestY = vertices[v];
					if (vertices[v] > largestY) largestY = vertices[v];
				} else {
					if (vertices[v] < smallestZ) smallestZ = vertices[v];
					if (vertices[v] > largestZ) largestZ = vertices[v];
				}
			}
//			std::cout << "Found smallest and largest" << std::endl;

			// move things so that the smallest vertex in all dim is at 1
			// scale everything down so that it fits from 0 to 1
			// might want to scale down even more to avoid the boundary situation
			T widthX = largestX - smallestX;
			T widthY = largestY - smallestY;
			T widthZ = largestZ - smallestZ;
			T scale = 1.f / (maxValueFromThree(widthX, widthY, widthZ) * 2.f);

			mesh_t mesh = shapes[0].mesh;

			int numVertices = vertices.size() / 3;
			double meshObjectPositions[vertices.size()];
			for (int i = 0; i < vertices.size(); i++) {
				// if (i % 3 == 0) meshObjectPositions[i] = (vertices[i] - smallestX) * scale;
				// else if (i % 3 == 1) meshObjectPositions[i] = (vertices[i] - smallestY) * scale;
				// else meshObjectPositions[i] = (vertices[i] - smallestZ) * scale;
				meshObjectPositions[i] = vertices[i];
			}
//			std::cout << "Moved and scale vertices" << std::endl;

			int numTriangles = mesh.indices.size() / 3;
			int meshObjectTriangles[mesh.indices.size()];
			for (int i = 0; i < mesh.indices.size(); i++) meshObjectTriangles[i] = mesh.indices[i].vertex_index;
//			std::cout << "Indices copied" << std::endl;

			// std::cout << "numVertices = " << numVertices << std::endl;
			// std::cout << "meshObjectPositions.size() = " << sizeof(meshObjectPositions) / sizeof(meshObjectPositions[0]) << std::endl;
			// std::cout << "numTriangles = " << numTriangles << std::endl;
			// std::cout << "meshObjectTriangles.size() = " << sizeof(meshObjectTriangles) / sizeof(meshObjectTriangles[0]) << std::endl;

			MeshObject* pMeshObject = construct_mesh_object(numVertices, &meshObjectPositions[0], numTriangles, &meshObjectTriangles[0]);
			std::cout << "Mesh object created" << std::endl;

			T estDistBwPars1D = dx / NpPerCell1D;
			for (int cx = offsetGrid; cx < res - 1 - offsetGrid; cx++) {
				for (int cy = offsetGrid; cy < res - 1 - offsetGrid; cy++) {
					for (int cz = offsetGrid; cz < res - 1 - offsetGrid; cz++) {
						TV cellWorldSpace = gridSpace2WorldSpace(cx, cy, cz);
						for (int px = 0; px < NpPerCell1D; px++) {
							for (int py = 0; py < NpPerCell1D; py++) {
								for (int pz = 0; pz < NpPerCell1D; pz++) {
									T offsetX = randWithOffset(0, estDistBwPars1D, 0.000001f);
									T offsetY = randWithOffset(0, estDistBwPars1D, 0.000001f);
									T offsetZ = randWithOffset(0, estDistBwPars1D, 0.000001f);

									TV parPos = TV::Zero();
									parPos[0] = cellWorldSpace[0] + px * estDistBwPars1D;// + offsetX;
									parPos[1] = cellWorldSpace[1] + py * estDistBwPars1D;// + offsetY;
									parPos[2] = cellWorldSpace[2] + pz * estDistBwPars1D;// + offsetZ;

									double testPoint[3] = {parPos[0], parPos[1], parPos[2]};
									if (point_inside_mesh(&testPoint[0], pMeshObject)) {
										xp.push_back(parPos);
									}
								}
							}
						}
					}
				}
			}
			destroy_mesh_object(pMeshObject);
//			std::cout << "Done sampling mesh" << std::endl;
		}
	}

	T maxValueFromThree(T a, T b, T c) {
		T maxVal = a;
		if (b <= a && c <= a) maxVal = a;
		if (a <= b && c <= b) maxVal = b;
		if (a <= c && b <= c) maxVal = c;
		return maxVal;
	}

	T randWithOffset(T lo, T hi, T offset) {
		T loWithOffset = lo + offset;
		T hiWithOffset = hi - offset;
		return ((hiWithOffset - loWithOffset) * ((T)rand() / RAND_MAX)) + loWithOffset;
	}

	TV gridSpace2WorldSpace(int x, int y, int z) {
		TV worldSpace = TV::Zero();
		worldSpace[0] = dx * x;
		worldSpace[1] = dx * y;
		worldSpace[2] = dx * z;
		return worldSpace;
	}

	int gridSpace2Idx(int x, int y, int z) {
		return x + y * res + z * res * res;
	}

	TV computeParticleMomentum(const std::vector<T> &mp, const std::vector<TV> &vp) {
		TV result = TV::Zero();
		for (int p = 0; p < Np; p++) {
			result += mp[p] * vp[p];
		}
		return result;
	}

	TV computeGridMomentum(const std::vector<T> &mg, const std::vector<TV>& vg) {
		TV result = TV::Zero();
		for (int i = 0; i < numNode; i++) {
			result += mg[i] * vg[i];
		}
		return result;
	}

	// Compute 1D quadratic B spline weights
	// x is assumed to be scaled in the index space (i.e., it is in a dx=1 grid)

	void computeWeights1D(TV& w, TV& dw, T& base_node, const T& x) {
		base_node = floor(x - 0.5) + 1;

		T d0 = x - base_node + 1;
		T z = 1.5 - d0;
		T z2 = z * z;
		w[0] = 0.5 * z2;

		T d1 = d0 - 1;
		w[1] = 0.75 - d1 * d1;

		T d2 = 1 - d1;
		T zz = 1.5 - d2;
		T zz2 = zz * zz;
		w[2] = 0.5 * zz2;

		dw[0] = -z;
		dw[1] = -2 * d1;
		dw[2] = zz;
	}

	void transferP2G(std::vector<T>& mg, std::vector<TV>& vg, std::vector<int>& active_nodes, const std::vector<TV>& xp, const std::vector<T>& mp, const std::vector<TV>& vp) {
		for (int p = 0; p < Np; p++) {
			TV X = xp[p];
			TV X_index_space = X / dx;
			TV w1, w2, w3, dw1, dw2, dw3;
			T base_node1, base_node2, base_node3;
			computeWeights1D(w1, dw1, base_node1, X_index_space[0]);
			computeWeights1D(w2, dw2, base_node2, X_index_space[1]);
			computeWeights1D(w3, dw3, base_node3, X_index_space[2]);

			if (base_node1 >= res|| base_node2 >= res|| base_node3 >= res || base_node1 < 0 || base_node2 < 0 || base_node3 < 0) {
				std::cout << "P2G: out of bound base nodes" << std::endl;
				std::cout << "base_node1 " << base_node1 << "; base_node2: " << base_node2 << "; base_node3: " << base_node3 << std::endl;
				std::cout << "X_index_space " << X_index_space[0] << " " << X_index_space[1] << " " << X_index_space[2] << std::endl;
				std::cout << "X " << X[0] << " " << X[1] << " " << X[2] << std::endl;
				std::getchar();
			}

			T sum_wijk = 0;
			for (int i = 0; i < 3; i++) {
				T wi = w1[i];
				T node_i = base_node1 + i;

				for (int j = 0; j < 3; j++) {
					T wij = wi * w2[j];
					T node_j = base_node2 + j;

					for (int k = 0; k < 3; k++) {
						T wijk = wij * w3[k];
						T node_k = base_node3 + k;

						// splat mass
						int idx = gridSpace2Idx(node_i, node_j, node_k);
						if (idx >= numNode || idx < 0) {
							std::cout << numNode << std::endl;
							std::cout << "P2G: idx out of bounds!" << std::endl;
							std::cout << "P2G: node_i = " << node_i << std::endl;
							std::cout << "P2G: node_j = " << node_j << std::endl;
							std::cout << "P2G: node_k = " << node_k << std::endl;
							std::cout << idx << std::endl;
							std::getchar();
						}
						mg[idx] += mp[p] * wijk;

						// splat momentum
						vg[idx] += wijk * mp[p] * vp[p];

						sum_wijk += wijk;
					}
				}
			}
		}

		for (int i = 0; i < numNode; i++) {
			if (mg[i] != 0) {
				active_nodes.push_back(i);
				vg[i] /= mg[i];
			} else {
				vg[i] = TV::Zero();
			}
		}
	}

	void addGravity(std::vector<TV>& force, const std::vector<T>& mg, const std::vector<int>& active_nodes, const TV& gravity) {
		for (int i = 0; i < active_nodes.size(); i++) {
			int idx = active_nodes[i];
			force[idx] += mg[idx] * gravity;
		}
	}

	void polarSVD(TSM& u, TSM& sigma, TSM& v, const TSM& F) {
		Eigen::JacobiSVD<TSM> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
		u = svd.matrixU();
		v = svd.matrixV();
		sigma = svd.singularValues().asDiagonal();

		if (u.determinant() < 0) {
			u(0, 2) *= -1;
			u(1, 2) *= -1;
			u(2, 2) *= -1;
			sigma(2, 2) *= -1;
		}

		if (v.determinant() < 0) {
			v(0, 2) *= -1;
			v(1, 2) *= -1;
			v(2, 2) *= -1;
			sigma(2, 2) *= -1;
		}
	}

	void fixedCorotated(const TSM& F, const T& mu, const T& lambda, TSM& P) {
		TSM u, sigma, v;
		polarSVD(u, sigma, v, F);
		TSM R = u * v.transpose();
		T J = F.determinant();
		TSM A = J * F.inverse().transpose();
		P = T(2) * mu * (F - R) + lambda * (J - 1) * A;
	}

	void addElasticity( std::vector<TV>& force, const std::vector<TV>& xp, const std::vector<TSM>& Fp, const T& Vp0, const T& mu, const T& lambda) {
		for (int p = 0; p < Np; p++) {
			TSM thisFp = Fp[p];
			TSM thisP;
			fixedCorotated(thisFp, mu, lambda, thisP);
			TSM Vp0PFt = Vp0 * thisP * thisFp.transpose();

			TV X = xp[p];
			TV X_index_space = X / dx;

			TV w1, w2, w3, dw1, dw2, dw3;
			T base_node1, base_node2, base_node3;

			computeWeights1D(w1, dw1, base_node1, X_index_space[0]);
			computeWeights1D(w2, dw2, base_node2, X_index_space[1]);
			computeWeights1D(w3, dw3, base_node3, X_index_space[2]);

			if (base_node1 >= res|| base_node2 >= res|| base_node3 >= res || base_node1 < 0 || base_node2 < 0 || base_node3 < 0) {
				std::cout << "addElasticity: out of bound base nodes" << std::endl;
				std::cout << "base_node1 " << base_node1 << "; base_node2: " << base_node2 << "; base_node3: " << base_node3 << std::endl;
				std::cout << "X_index_space: " << X_index_space[0] << " " << X_index_space[1] << " " << X_index_space[2] << std::endl;
				std::getchar();
			}

			for (int i = 0; i < 3; i++) {
				T wi = w1[i];
				T dwidxi = dw1[i] / dx;
				T node_i = base_node1 + i;

				for (int j = 0; j < 3; j++) {
					T wj = w2[j];
					T wij = wi * wj;
					T dwijdxi = dwidxi * wj;
					T dwijdxj = wi / dx * dw2[j];
					T node_j = base_node2 + j;

					for (int k = 0; k < 3; k++) {
						T wk = w3[k];
						T wijk = wij * wk;
						T dwijkdxi = dwijdxi * wk;
						T dwijkdxj = dwijdxj * wk;
						T dwijkdxk = wij / dx * dw3[k];
						T node_k = base_node3 + k;

						TV grad_w = TV(dwijkdxi, dwijkdxj, dwijkdxk);
						TV foo = -Vp0PFt * grad_w;

						int idx = gridSpace2Idx(node_i, node_j, node_k);
						if (idx >= numNode || idx < 0) {
							std::cout << "addElasticity: idx out of bound!" << std::endl;
							std::cout << "node_i = " << node_i << std::endl;
							std::cout << "node_j = " << node_j << std::endl;
							std::cout << "node_k = " << node_k << std::endl;
							std::getchar();
						}
						force[idx] += foo;
					}
				}
			}
		}
	}

	void updateGridVelocity(const std::vector<T>& mg, const std::vector<TV>& vgn, const std::vector<TV>& force, const std::vector<int>& active_nodes, const T& dt, std::vector<TV>& vg) {
		for (int i = 0; i < active_nodes.size(); i++) {
			int idx = active_nodes[i];
			vg[idx] = vgn[idx] + dt * force[idx] / mg[idx];
		}
	}

	void setBoundaryVelocities(int thickness, std::vector<TV>& vg) {
		// x direction
		for (int i = 0; i < thickness; i ++) {
			for (int j = 0; j < res; j++) {
				for (int k = 0; k < res; k++) {
					vg[gridSpace2Idx(i, j, k)] = TV::Zero();
				}
			}
		}

		for (int i = res - thickness; i < res; i ++) {
			for (int j = 0; j < res; j++) {
				for (int k = 0; k < res; k++) {
					vg[gridSpace2Idx(i, j, k)] = TV::Zero();
				}
			}
		}

		// y direction
		for (int i = 0; i < res; i ++) {
			for (int j = 0; j < thickness; j++) {
				for (int k = 0; k < res; k++) {
					vg[gridSpace2Idx(i, j, k)] = TV::Zero();
				}
			}
		}

		for (int i = 0; i < res; i ++) {
			for (int j = res - thickness; j < res; j++) {
				for (int k = 0; k < res; k++) {
					vg[gridSpace2Idx(i, j, k)] = TV::Zero();
				}
			}
		}

		// z direction
		for (int i = 0; i < res; i ++) {
			for (int j = 0; j < res; j++) {
				for (int k = 0; k < thickness; k++) {
					vg[gridSpace2Idx(i, j, k)] = TV::Zero();
				}
			}
		}

		for (int i = 0; i < res; i ++) {
			for (int j = 0; j < res; j++) {
				for (int k = res - thickness; k < res; k++) {
					vg[gridSpace2Idx(i, j, k)] = TV::Zero();
				}
			}
		}
	}

	void transferG2P(const T& dt, const std::vector<TV>& vgn, const std::vector<TV>& vg, const T& flip, std::vector<TV>& xp, std::vector<TV>& vp) {
		for (int p = 0; p < Np; p++) {
			TV X = xp[p];
			TV X_index_space = X / dx;
			TV w1, w2, w3, dw1, dw2, dw3;
			T base_node1, base_node2, base_node3;

			computeWeights1D(w1, dw1, base_node1, X_index_space[0]);
			computeWeights1D(w2, dw2, base_node2, X_index_space[1]);
			computeWeights1D(w3, dw3, base_node3, X_index_space[2]);

			if (base_node1 >= res|| base_node2 >= res|| base_node3 >= res || base_node1 < 0 || base_node2 < 0 || base_node3 < 0) {
				std::cout << "G2P: out of bound base nodes" << std::endl;
				std::cout << "base_node1 " << base_node1 << "; base_node2: " << base_node2 << "; base_node3: " << base_node3 << std::endl;
				std::cout << "X_index_space: " << X_index_space[0] << " " << X_index_space[1] << " " << X_index_space[2] << std::endl;
				std::getchar();
			}

			TV vpic = TV::Zero();
			TV vflip = vp[p];

			T sum_wijk = 0;

			for (int i = 0; i < 3; i++) {
				T wi = w1[i];
				T node_i = base_node1 + i;

				for(int j = 0; j < 3; j++) {
					T wij = wi * w2[j];
					T node_j = base_node2 + j;

					for (int k = 0; k < 3; k++) {
						T wijk = wij * w3[k];
						T node_k = base_node3 + k;

						int idx = gridSpace2Idx(node_i, node_j, node_k);
						if (idx >= numNode || idx < 0) {
							std::cout << "G2P: idx out of bound!" << std::endl;
							std::cout << "node_i = " << node_i << std::endl;
							std::cout << "node_j = " << node_j << std::endl;
							std::cout << "node_k = " << node_k << std::endl;
							std::getchar();
						}
						vpic += wijk * vg[idx];
						vflip += wijk * (vg[idx] - vgn[idx]);
						sum_wijk += wijk;
					}
				}
			}
			vp[p] = (1 - flip) * vpic + flip * vflip;
			xp[p] += dt * vpic;
		}
	}

	void evolveF(const T& dt, const std::vector<TV>& vg, const std::vector<TV>& xp, std::vector<TSM>& Fp) {
		for (int p = 0; p < Np; p++) {
			TSM thisFp = Fp[p];

			TV X = xp[p];
			TV X_index_space = X / dx;

			TV w1, w2, w3, dw1, dw2, dw3;
			T base_node1, base_node2, base_node3;

			computeWeights1D(w1, dw1, base_node1, X_index_space[0]);
			computeWeights1D(w2, dw2, base_node2, X_index_space[1]);
			computeWeights1D(w3, dw3, base_node3, X_index_space[2]);

			if (base_node1 >= res|| base_node2 >= res|| base_node3 >= res || base_node1 < 0 || base_node2 < 0 || base_node3 < 0) {
				std::cout << "evolvF: out of bound base nodes" << std::endl;
				std::cout << "base_node1 " << base_node1 << "; base_node2: " << base_node2 << "; base_node3: " << base_node3 << std::endl;
				std::cout << "X_index_space: " << X_index_space[0] << " " << X_index_space[1] << " " << X_index_space[2];
				std::getchar();
			}

			// Compute grad_vp
			TSM grad_vp = TSM::Zero(dim, dim);
			for (int i = 0; i < 3; i++) {
				T wi = w1[i];
				T dwidxi = dw1[i] / dx;
				T node_i = base_node1 + i;

				for (int j = 0; j < 3; j++) {
					T wj = w2[j];
					T wij = wi * wj;
					T dwijdxi = dwidxi * wj;
					T dwijdxj = wi / dx * dw2[j];
					T node_j = base_node2 + j;

					for (int k = 0; k < 3; k++) {
						T wk = w3[k];
						T wijk = wij * wk;
						T dwijkdxi = dwijdxi * wk;
						T dwijkdxj = dwijdxj * wk;
						T dwijkdxk = wij / dx * dw3[k];
						T node_k = base_node3 + k;

						TV grad_w = TV(dwijkdxi, dwijkdxj, dwijkdxk);
						int idx = gridSpace2Idx(node_i, node_j, node_k);
						if (idx >= numNode || idx < 0) {
							std::cout << "evolveF: idx out of bound!" << std::endl;
							std::cout << "node_i = " << node_i << std::endl;
							std::cout << "node_j = " << node_j << std::endl;
							std::cout << "node_k = " << node_k << std::endl;
							std::getchar();
						}
						TV vijk = vg[idx];
						grad_vp += vijk * grad_w.transpose();
					}
				}
			}

			TSM newFp = (TSM::Identity() + dt * grad_vp) * thisFp;

			for (int i = 0; i < dim; i++) {
				for (int j = 0; j < dim; j++) {
					Fp[p](i, j) = newFp(i, j);
				}
			}
		}
	}

	void advanceOneStep() {
		// Init zero grid data
		std::vector<T> mg = std::vector<T>(numNode, 0); // grid mass
		std::vector<TV> vgn = std::vector<TV>(numNode, TV::Zero()); // grid new velocity
		std::vector<TV> vg = std::vector<TV>(numNode, TV::Zero()); // grid old velocity
		std::vector<TV> force = std::vector<TV>(numNode, TV::Zero());
		std::vector<int> active_nodes;

		// P2G
		//TV Lp = computeParticleMomentum(mp, vp);
		transferP2G(mg, vgn, active_nodes, xp, mp, vp);

		// Compute force
		addGravity(force, mg, active_nodes, gravity);
		addElasticity(force, xp, Fp, Vp0, mu, lambda);

		// Update velocity
		updateGridVelocity(mg, vgn, force, active_nodes, dt, vg);

		// Boundary conditions
		setBoundaryVelocities(5, vg);

		// G2P
		// Lg = computeGridMomentum(mg, vg);
		evolveF(dt, vg, xp, Fp);
		transferG2P(dt, vgn, vg, 0.95, xp, vp);
		// Lp = computeParticleMomentum(mp, vp);
	}

	void dumpPoly(std::string filename)
	{
		std::ofstream fs;
		fs.open(filename);
		fs << "POINTS\n";
		int count = 0;
		for (auto parPos : xp) {
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

	void run(const int max_frame)
	{
		for(int frame=1; frame<max_frame; frame++) {
			//for(int frame=1; frame<5; frame++) {
			std::cout << "Frame " << frame << std::endl;

			int N_substeps = (int)(((T)1/24)/dt);
			for (int step = 1; step <= N_substeps; step++) {
				//for (int step = 1; step <= 5; step++) {
				std::cout << "Step " << step << std::endl;
				advanceOneStep();
			}
			mkdir("output/", 0777);
			std::string filename = "output/" + std::to_string(frame) + ".poly";
			dumpPoly(filename);
			std::cout << std::endl;
		}
	}
};
