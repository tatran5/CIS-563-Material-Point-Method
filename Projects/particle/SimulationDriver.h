
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

	ParticleSystem<T,dim> mPs;
	T mDt;
	TV mGravity;
	T mGround;
	T mFlipWeight;

	SimulationDriver()
    // : dt((T)0.00001)
    : mDt((T)0.001),  // 150 times bigger dt than explicit. We can't go arbitrarily large because we are still doing approximations to the non-linear problem using taylor expansion.
    mFlipWeight((T)0.95)
    {
    	mGravity.setZero();
    	mGravity(1) = -9.8;
    	mGround = 0.1;
    }

    void run(const int max_frame)
    {
    	for(int frame=1; frame<max_frame; frame++) {
    	//for(int frame=1; frame<5; frame++) {
    		std::cout << "Frame " << frame << std::endl;

    		int N_substeps = (int)(((T)1/24)/mDt);
    		for (int step = 1; step <= N_substeps; step++) {
    	 	//for (int step = 1; step <= 5; step++) {
    			std::cout << "Step " << step << std::endl;
    			advanceOneStep();
    			//std::cout << std::endl;
    		//	mkdir("output/", 0777);
    		// std::string filename = "output/a_" + std::to_string(step) + ".poly";
    		// mPs.dumpPoly(filename);
    		// std::cout << std::endl;
    		}
    		mkdir("output/", 0777);
    		std::string filename = "output/" + std::to_string(frame) + ".poly";
    		mPs.dumpPoly(filename);
    		std::cout << std::endl;
    	}
    }

    void advanceOneStep() {
    // Clean grid data: mass, position, velocity
    	int numNode = mPs.mNodeVel.size();
    	mPs.mNodeMass = std::vector<T>(numNode, 0.f);
    	// mPs.mNodePos = std::vector<TV>(numNode, TV::Zero());
    	mPs.mNodeVel = std::vector<TV>(numNode, TV::Zero());
    	std::vector<TV> nodeVelNew(numNode, TV::Zero());
    	std::vector<TV> forceOnNode(numNode, TV::Zero());
    	std::vector<TV> activeNodes = std::vector<TV>();

    	// P2G -----
    	//TV totalParMomentumP2G = mPs.getTotalParticleMomentum();
    	mPs.transferP2G(&activeNodes, &nodeVelNew);
    	//TV totalGridMomentumP2G = mPs.getTotalGridMomentum(mPs.mNodeMass, nodeVelNew);
    	// std::cout << "totalParMomentumP2G " << std::endl << totalParMomentumP2G << std::endl;
    	// std::cout << "totalGridMomentumP2G " << std::endl << totalGridMomentumP2G << std::endl;

 	    // Add force -----
    	mPs.computeGravityOnGrid(&forceOnNode, activeNodes, mGravity);
    	mPs.computeElasticityOnGrid(&forceOnNode);

    	// Update grid velocity -----
    	mPs.updateGridVelocity(activeNodes, nodeVelNew, forceOnNode, mGravity, mDt);

    	// Set boundary conditions
    	int thickness = 1;
    	mPs.setBoundaryVelocities(thickness);

    	// Transfer G2P -----
		// TV totalGridMomentumG2P = mPs.getTotalGridMomentum(mPs.mNodeMass, mPs.mNodeVel);
		mPs.evolveF(mDt);
		mPs.transferG2P(nodeVelNew, mFlipWeight, mDt);
      	// TV totalParMomentumG2P = mPs.getTotalParticleMomentum();
		// std::cout << "totalGridMomentumG2P " << std::endl << totalGridMomentumG2P << std::endl;
		// std::cout << "totalParMomentumG2P " << std::endl << totalParMomentumG2P << std::endl; 
    }
};
