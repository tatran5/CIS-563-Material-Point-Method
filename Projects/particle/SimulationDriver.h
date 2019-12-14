#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "ParticleSystem.h"

using namespace std;

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
	vector<TV> xp; // position of particle
	vector<T> mp; // mass of particle
	vector<TV> vp; // velocity of particle
	vector<TSM> Fp;

	// Other
	TSM identity;

	SimulationDriver() {
		// Set values for simulation parameters
		dt = 1e-3;
		gravity = TV(0, -9.8, 0);

		// Set values for grid parameters
		// dx = 0.02f;
		dx = 0.05f;
		minCorner = 0.f;
		maxCorner = 1.f;
		res = (maxCorner - minCorner) / dx + 1;
		// Set values for add-on grid parameters
		NpPerCell1D = 2.f;
		offsetGrid = 4.f;
		numNode = res * res * res;

		// Set values for particles
		E = 10000;
		nu = 0.3f;
		mu = E / (2 * (1 + nu));
		lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
		rho = 1000;
		sampleParticles();
		Np = xp.size();
		Vp0 = dx * dx * dx / NpPerCell1D / NpPerCell1D / NpPerCell1D;
		mp = vector<T>(Np, Vp0 * rho);
		vp = vector<TV>(Np, TV::Zero());
		identity.setIdentity();
		Fp = vector<TSM>(Np, identity);
	}

	void sampleParticles() {
		T estDistBwPars1D = dx / NpPerCell1D;
		for (int cx = offsetGrid; cx < res - 1 - offsetGrid; cx++) {
			for (int cy = offsetGrid; cy < res - 1 - offsetGrid; cy++) {
				for (int cz = offsetGrid; cz < res - 1 - offsetGrid; cz++) {
					TV cellWorldSpace = gridSpace2WorldSpace(cx, cy, cz);
					// cout << "c0 " << cellWorldSpace[0] << endl;
					// cout << "c1 " << cellWorldSpace[1] << endl;
					// cout << "c2 " << cellWorldSpace[2] << endl;
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
								// cout << "p0 " << parPos[0] << endl;
								// cout << "p1 " << parPos[1] << endl;
								// cout << "p2 " << parPos[2] << endl;
							}
						}
					}
				}
			}
		}
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

	TV computeParticleMomentum(const vector<T> &mp, const vector<TV> &vp) {
		TV result = TV::Zero();
		for (int p = 0; p < Np; p++) {
			result += mp[p] * vp[p];
		}
		return result;
	}

	TV computeGridMomentum(const vector<T> &mg, const vector<TV>& vg) {
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

	void transferP2G(vector<T>& mg, vector<TV>& vg, vector<int>& active_nodes, const vector<TV>& xp, const vector<T>& mp, const vector<TV>& vp) {
		for (int p = 0; p < Np; p++) {
			TV X = xp[p];
			TV X_index_space = X / dx;
			// cout << X_index_space << endl;
			TV w1, w2, w3, dw1, dw2, dw3;
			T base_node1, base_node2, base_node3;
			computeWeights1D(w1, dw1, base_node1, X_index_space[0]);
			computeWeights1D(w2, dw2, base_node2, X_index_space[1]);
			computeWeights1D(w3, dw3, base_node3, X_index_space[2]);

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
						if (idx > numNode) {
							cout << "P2G: idx out of bounds!" << endl;
							cout << "P2G: node_i = " << node_i << endl;
							cout << "P2G: node_j = " << node_j << endl;
							cout << "P2G: node_k = " << node_k << endl;
						}
						mg[idx] += mp[p] * wijk;

						// splat momentum
						vg[idx] += wijk * mp[p] * vp[p];

						sum_wijk += wijk;
					}
				}
			}
			// cout << "sum_wijk = " << sum_wijk << endl;
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

	void addGravity(vector<TV>& force, const vector<T>& mg, const vector<int>& active_nodes, const TV& gravity) {
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
		TSM FChanged = u * sigma * v.transpose();
		TSM R = u * v.transpose();
		T J = FChanged.determinant();
		//cout << "J = " << J << endl;

		TSM A = FChanged.adjoint().transpose();
		P = 2 * mu * (F - R) + lambda * (J - 1) * A;
	}

	void addElasticity( vector<TV>& force, const vector<TV>& xp, const vector<TSM>& Fp, const T& Vp0, const T& mu, const T& lambda) {
		//vector<TV> forceBefore = force;
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
						if (idx > numNode) {
							cout << "G2P: idx out of bound!" << endl;
							cout << "node_i = " << node_i << endl;
							cout << "node_j = " << node_j << endl;
							cout << "node_k = " << node_k << endl;
						}
						// cout << "foo = " << endl << foo << endl;
						// cout << "grad_w = " << endl << grad_w << endl;
						// cout << "Vp0PFt = " << endl << Vp0PFt << endl;
						// cout << "Vp0 = " << endl << Vp0 << endl;
						force[idx] += foo;
					}
				}
			}
		}
	}

	void updateGridVelocity(const vector<T>& mg, const vector<TV>& vgn, const vector<TV>& force, const vector<int>& active_nodes, const T& dt, vector<TV>& vg) {
		for (int i = 0; i < active_nodes.size(); i++) {
			int idx = active_nodes[i];
			vg[idx] = vgn[idx] + dt * force[idx] / mg[idx];
		}
	}

	void setBoundaryVelocities(int thickness, vector<TV>& vg) {
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

	void transferG2P(const T& dt, const vector<TV>& vgn, const vector<TV>& vg, const T& flip, vector<TV>& xp, vector<TV>& vp) {
		for (int p = 0; p < Np; p++) {
			TV X = xp[p];
			TV X_index_space = X / dx;
			TV w1, w2, w3, dw1, dw2, dw3;
			T base_node1, base_node2, base_node3;

			computeWeights1D(w1, dw1, base_node1, X_index_space[0]);
			computeWeights1D(w2, dw2, base_node2, X_index_space[1]);
			computeWeights1D(w3, dw3, base_node3, X_index_space[2]);

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
						if (idx > numNode) {
							cout << "G2P: idx out of bound!" << endl;
							cout << "node_i = " << node_i << endl;
							cout << "node_j = " << node_j << endl;
							cout << "node_k = " << node_k << endl;
						}
						vpic += wijk * vg[idx];
						vflip += wijk * (vg[idx] - vgn[idx]);
						sum_wijk += wijk;
					}
				}
			}
			// cout << "sum_wijk = " << sum_wijk << endl;
			vp[p] = (1 - flip) * vpic + flip * vflip;
			xp[p] += dt * vpic;
		}
	}

	void evolveF(const T& dt, const vector<TV>& vg, const vector<TV>& xp, vector<TSM>& Fp) {
		for (int p = 0; p < Np; p++) {
			TSM thisFp = Fp[p];

			TV X = xp[p];
			TV X_index_space = X / dx;

			TV w1, w2, w3, dw1, dw2, dw3;
			T base_node1, base_node2, base_node3;

			computeWeights1D(w1, dw1, base_node1, X_index_space[0]);
			computeWeights1D(w2, dw2, base_node2, X_index_space[1]);
			computeWeights1D(w3, dw3, base_node3, X_index_space[2]);

			// Compute grad_vp
			TSM grad_vp;
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
						if (idx > numNode) {
							cout << "G2P: idx out of bound!" << endl;
							cout << "node_i = " << node_i << endl;
							cout << "node_j = " << node_j << endl;
							cout << "node_k = " << node_k << endl;
						}
						TV vijk = vg[idx];
						grad_vp += vijk * grad_w.transpose();
					}
				}
			}

			TSM newFp = (identity + dt * grad_vp) * thisFp;

			for (int i = 0; i < dim; i++) {
				for (int j = 0; j < dim; j++) {
					Fp[p](i, j) = newFp(i, j);
				}
			}
		}
	}

	void advanceOneStep() {
		// Init zero grid data
		vector<T> mg = vector<T>(numNode, 0); // grid mass
		vector<TV> vgn = vector<TV>(numNode, TV::Zero()); // grid new velocity
		vector<TV> vg = vector<TV>(numNode, TV::Zero()); // grid old velocity
		vector<TV> force = vector<TV>(numNode, TV::Zero());
		vector<int> active_nodes;

		// P2G
		//TV Lp = computeParticleMomentum(mp, vp);
		transferP2G(mg, vgn, active_nodes, xp, mp, vp);
		//TV Lg = computeGridMomentum(mg, vgn);
		// cout << "P2G: Lp = " << endl << Lp << endl;
		// cout << "P2G: Lg = " << endl << Lg << endl;

		// Compute force
		addGravity(force, mg, active_nodes, gravity);
		addElasticity(force, xp, Fp, Vp0, mu, lambda);

		// Update velocity
		updateGridVelocity(mg, vgn, force, active_nodes, dt, vg);
		//
		// Boundary conditions
		setBoundaryVelocities(1, vg);

		// G2P
		// Lg = computeGridMomentum(mg, vg);
		evolveF(dt, vg, xp, Fp);
		transferG2P(dt, vgn, vg, 0.95, xp, vp);
		// Lp = computeParticleMomentum(mp, vp);
		// cout << "G2P: Lp = " << endl << Lp << endl;
		// cout << "G2P: Lg = " << endl << Lg << endl;
	}

	void dumpPoly(string filename)
	{
		ofstream fs;
		fs.open(filename);
		fs << "POINTS\n";
		int count = 0;
		for (auto parPos : xp) {
			fs << ++count << ":";
			for (int i = 0; i < dim; i++){
				if (parPos(i) != parPos(i))
				cout << "GETTING NAN FOR PARTICLE POSITION" << endl;
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
			cout << "Frame " << frame << endl;

			int N_substeps = (int)(((T)1/24)/dt);
			for (int step = 1; step <= N_substeps; step++) {
				//for (int step = 1; step <= 5; step++) {
				// cout << "Step " << step << endl;
				advanceOneStep();
			}
			mkdir("output/", 0777);
			string filename = "output/" + to_string(frame) + ".poly";
			dumpPoly(filename);
			cout << endl;
		}
	}
};
