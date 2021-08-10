#pragma once
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include "ichnos_structures.h"

namespace ICHNOS{

    class XYZ_base{
    public:
        XYZ_base(boost::mpi::communicator& world_in);
        virtual void readXYZdata(std::string xyz_file){}
        virtual void calcWeights(vec3& p, std::vector<int>& ids, std::vector<double>& weights){}
        virtual void reset(){}

    protected:
        boost::mpi::communicator world;
    };

    XYZ_base::XYZ_base(boost::mpi::communicator& world_in)
        :
        world(world_in)
    {}
}
