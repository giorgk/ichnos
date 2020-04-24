#pragma once
#include "ichnos_options.h"
#include "ichnos_structures.h"
#include "velocity_base.h"


namespace ICHNOS {

	class layervelocity {
	public:
		layervelocity() {};
		std::vector<double> elevation;
		std::vector< std::vector<double> > vx;
		std::vector< std::vector<double> > vy;
		std::vector< std::vector<double> > vz;
	};

	class layer3D : public velocityField {
	public:
		layer3D(boost::mpi::communicator& world_in);
		void readVelocityField(std::string vf_file);
		void calcVelocity(vec3& vel, std::map<int, double>& proc_map, vec3& p);

	private:
		// point cloud containing elevation and velocity defined on the same points
		pointCloud< layervelocity > pcEV;

	};

	layer3D::layer3D(boost::mpi::communicator& world_in)
		:
		velocityField(world_in)
	{}


	void layer3D::readVelocityField(std::string vf_file) {

	}

	void layer3D::calcVelocity(vec3& vel, std::map<int, double>& proc_map, vec3& p) {

	}
}
