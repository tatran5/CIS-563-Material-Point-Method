
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "ParticleSystem.h"

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

	// Particle parameters
	T E; //elasticity
	T nu, mu, lambda, rho, Np, Vp0;
	std::vector<TV> xp; // position of particle
	std::vector<T> mp; // mass of particle
	std::vector<TV> vp; // velocity of particle
	std::vector<TSM> Fp;

	SimulationDriver() {
		// Set values for simulation parameters
		dt = 1e-3;
		gravity = TV(0, -9.8, 0);

		// Set values for grid parameters
		// dx = 0.02f;
		dx = 0.1f;
		minCorner = 0.f;
		maxCorner = 1.f;
		res = (maxCorner - minCorner) / dx + 1;
		// Set values for add-on grid parameters
		NpPerCell1D = 2.f;
		offsetGrid = 3.f;

		// Set values for particles
		E = 1e4;
		nu = 0.3f;
		mu = E / (2 * (1 + nu));
		lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
		rho = 1000;
		Np = xp.size();
		Vp0 = dx * dx * dx / NpPerCell1D / NpPerCell1D / NpPerCell1D;
		mp = std::vector<T>(Np, Vp0 * rho);
		vp = std::vector<TV>(Np, TV::Zero());
		TSM identity;
		identity.setIdentity();
		Fp = std::vector<TSM>(Np, identity);

		sampleParticles();
	}

	void sampleParticles() {
		T estDistBwPars1D = dx / NpPerCell1D;
		for (int cx = offsetGrid; cx < res - 1 - offsetGrid; cx++) {
			for (int cy = offsetGrid; cy < res - 1 - offsetGrid; cy++) {
				for (int cz = offsetGrid; cz < res - 1 - offsetGrid; cz++) {
					TV cellWorldSpace = gridSpace2WorldSpace(cx, cy, cz);
					// std::cout << "c0 " << cellWorldSpace[0] << std::endl;
					// std::cout << "c1 " << cellWorldSpace[1] << std::endl;
					// std::cout << "c2 " << cellWorldSpace[2] << std::endl;
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
								// std::cout << "p0 " << parPos[0] << std::endl;
								// std::cout << "p1 " << parPos[1] << std::endl;
								// std::cout << "p2 " << parPos[2] << std::endl;
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
		return ((hiWithOffset - loWithOffset) * ((T)std::rand() / RAND_MAX)) + loWithOffset;
	}


	TV gridSpace2WorldSpace(int x, int y, int z) {
		TV worldSpace = TV::Zero();
		worldSpace[0] = dx * x;
		worldSpace[1] = dx * y;
		worldSpace[2] = dx * z;
		return worldSpace;
	}

	TV computeParticleMomentum(const std::vector<T> &mp, const std::vector<TV> &vp) {
		TV result = TV::Zero();
		for (int p = 0; p < Np; p++) {
			result += mp[p] * vp[p];
		}
		return result;
	}

	void advanceOneStep() {
		// Init zero grid data
		int numNode = res * res * res;
		std::vector<T> mg = std::vector<T>(numNode, 0); // grid mass
		std::vector<TV> vgn = std::vector<TV>(numNode, TV::Zero()); // grid new velocity
		std::vector<TV> vg = std::vector<TV>(numNode, TV::Zero()); // grid old velocity
		std::vector<TV> force = std::vector<TV>(numNode, TV::Zero());

		// P2G
		TV Lp = computeParticleMomentum(mp, vp);
		transferP2G(vgn);

	}

	// X is location of a particle in grid space
	void computeWeights1D(TV& w, T& base_node, const T& x) {
		T base_node = std::floor(x - 0.5f) + 1;

	}

	void transferP2G(std::vector<TV>& vg) {
		for (int p = 0; p < Np; p++) {
			TV X = xp[p];
			TV X_index_space = X / dx;
			TV w1, w2, w3;
			T base_node1, base_node2, base_node3;
			computeWeights1D(X_index_space[0]);
			computeWeights1D(X_index_space[1]);
			computeWeights1D(X_index_space[2]);

			for (int i = 0; i < 3; i++) {
				T wi = w1[i];
				T node_i = base_node1 + i;

				for (int j = 0; j < 3; j++) {
					T wij = wi * w2[j];
					T node_j = base_node2 + j;

					for (int k = 0; k < 3; k++) {
						T wijk = wij * w3[k];

					}
				}
			}
		}
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
		//for(int frame=1; frame<max_frame; frame++) {
		for(int frame=1; frame<5; frame++) {
			std::cout << "Frame " << frame << std::endl;

			int N_substeps = (int)(((T)1/24)/dt);
			//for (int step = 1; step <= N_substeps; step++) {
			for (int step = 1; step <= 5; step++) {
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
