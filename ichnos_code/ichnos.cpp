// ichnos.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <boost/mpi.hpp>

#include <iostream>

#include "ichnos_options.h"
#include "velocity_base.h"
#include "trace.h"

int main(int argc, char* argv[])
{
    boost::mpi::environment env{ argc, argv };
    boost::mpi::communicator world;

    if (world.rank() == 0) {
        std::cout << "Ichnos will run using " << world.size() << " processors" << std::endl;
        std::cout << "Good luck with that!" << std::endl;

    }

    ICHNOS::options OPT(world);
    if (!OPT.readInput(argc, argv))
        return 0;

    ICHNOS::pc3D VF(world);
    VF.readVelocityField(OPT.getVelFname());

    ICHNOS::Domain2D domain(OPT.Dopt);
    
    std::cout << world.rank() << " has " << VF.getCloudSize() << " points" << std::endl;

    ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);

    pt.Trace();

    return 0;
 }
