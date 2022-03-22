//
// Created by giorg on 3/21/2022.
//

#ifndef ICHNOS_ICHNOS_RWPT_H
#define ICHNOS_ICHNOS_RWPT_H

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include "ichnos_structures.h"
#include "TransientVelocity.h"

namespace ic = ICHNOS;
namespace ICHNOS{
    class CloudRWVel : public ic::CloudVel{
    public:
        CloudRWVel(boost::mpi::communicator& world_in, XYZ_base &XYZ_in);
        bool readVelocityField(std::string vf_file, int nPnts);
        void calcVelocity(vec3& p, vec3& vel,
                          std::map<int, double>& proc_map,
                          helpVars& pvlu,
                          bool& out,
                          double time = 0);
        //void reset();
        //void updateStep(double& step);
        //void getVec3Data(std::vector<ic::vec3>& data);

    private:
        double aL_Scalar;
        double aT_Scalar;
        double Dm_Scalar;
        double Dt;
        double Ds;
        ICHNOS::SingletonGenerator* RG = RG->getInstance();

    };

    CloudRWVel::CloudRWVel(boost::mpi::communicator &world_in, XYZ_base &XYZ_in)
        :
            CloudVel(world_in,XYZ_in)
    {
        InterpolateOutsideDomain = true;
    }

    bool CloudRWVel::readVelocityField(std::string vf_file, int nPnts) {
        CloudVel::readVelocityField(vf_file, nPnts);
        if (world.rank() == 0)
            std::cout << "--> Random Walk related parameters file: " << vf_file << std::endl;

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
        ("RWPT.aL", po::value<std::string>(), "Longitudinal Dispersivity")
        ("RWPT.aT", po::value<std::string>(), "Transversal Dispersivity")
        ("RWPT.Dm", po::value<std::string>(), "Molecular diffusion")
        ("RWPT.Dt", po::value<double>()->default_value(1.0), "Delta Time")
        ("RWPT.Ds", po::value<double>()->default_value(1.0), "Derivative distance")
        ;

        po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions, true), vm_vfo);

        { //RWPT options
            std::string tmp_name;
            tmp_name = vm_vfo["RWPT.aL"].as<std::string>();
            if (is_input_scalar(tmp_name)){
                aL_Scalar = std::stof(tmp_name);
            }
            tmp_name = vm_vfo["RWPT.aT"].as<std::string>();
            if (is_input_scalar(tmp_name)){
                aT_Scalar = std::stof(tmp_name);
            }
            tmp_name = vm_vfo["RWPT.Dm"].as<std::string>();
            if (is_input_scalar(tmp_name)){
                Dm_Scalar = std::stof(tmp_name);
            }
            Dt = vm_vfo["RWPT.Dt"].as<double>();
            Ds = vm_vfo["RWPT.Ds"].as<double>();
        }
        return true;
    }



    void CloudRWVel::calcVelocity(vec3 &p, vec3 &vel,
                                  std::map<int, double> &proc_map,
                                  helpVars &pvlu,
                                  bool &out,
                                  double tm) {

        double aL = aL_Scalar;
        double aT = aT_Scalar;
        double Dm = Dm_Scalar;

        // Velocity along X
        vec3 px1 = p, vx1;
        vec3 px2 = p, vx2;
        px1.x = p.x + Ds;
        px2.x = p.x - Ds;
        CloudVel::calcVelocity(px1, vx1, proc_map, pvlu, out, tm);
        CloudVel::calcVelocity(px2, vx2, proc_map, pvlu, out, tm);

        // Velocity along Y
        vec3 py1 = p, vy1;
        vec3 py2 = p, vy2;
        py1.y = p.y + Ds;
        py2.y = p.y - Ds;
        CloudVel::calcVelocity(py1, vy1, proc_map, pvlu, out, tm);
        CloudVel::calcVelocity(py2, vy2, proc_map, pvlu, out, tm);

        // Velocity along Y
        vec3 pz1 = p, vz1;
        vec3 pz2 = p, vz2;
        pz1.z = p.z + Ds;
        pz2.z = p.z - Ds;
        CloudVel::calcVelocity(pz1, vz1, proc_map, pvlu, out, tm);
        CloudVel::calcVelocity(pz2, vz2, proc_map, pvlu, out, tm);

        vec3 v;
        CloudVel::calcVelocity(p, v, proc_map, pvlu, out, tm);

        // Calculate velocity derivatives
        double vx_partial = (vx1.x - vx2.x)/(2*Ds);
        double vy_partial = (vy1.y - vy2.y)/(2*Ds);
        double vz_partial = (vz1.z - vz2.z)/(2*Ds);

        double vlen = v.len();
        double v2xy = std::sqrt(v.x * v.x + v.y * v.y);
        double v2x = v.x*v.x;
        double v2y = v.y*v.y;
        double v2z = v.z*v.z;
        double v3 = vlen * vlen * vlen;
        // Calculate derivatives of dispersion coefficients
        double Dxx_p = v.x * vx_partial * (aL*(2/vlen - v2x / v3) - aT*((v2y+v2z)/v3));
        double Dxy_p = (aL - aT) * (vy_partial * v.x / vlen - vy_partial * v.x * v2y / v3);
        double Dxz_p = (aL - aT) * (vz_partial * v.x / vlen - vz_partial * v.x * v2z / v3);

        double Dyy_p = v.y * vy_partial * (aL * (2 / vlen - v2y / v3) - aT * ((v2x + v2z) / v3));
        double Dyx_p = (aL - aT) * (vx_partial * v.y / vlen - vx_partial * v.y * v2x / v3);
        double Dyz_p = (aL - aT) * (vz_partial * v.y / vlen - vz_partial * v.y * v2z / v3);

        double Dzz_p = v.z * vz_partial * (aL * (2 / vlen - (v2z) / v3) - aT * ((v2x + v2y) / v3));
        double Dzx_p = (aL - aT) * (vx_partial * v.z / vlen - vx_partial * v.z * v2x / v3);
        double Dzy_p = (aL - aT) * (vy_partial * v.z / vlen - vy_partial * v.z * v2y / v3);

        //double DL = aL*vlen;
        //double DT = aT*vlen;
        double sqaL = std::sqrt(2*(aL*vlen + Dm));
        double sqaT = std::sqrt(2*(aT*vlen + Dm));
        double sqDTime = std::sqrt(Dt);
/*
        double Dxx = aT*vlen + (aL - aT)*(v2x/vlen) + Dm;
        double Dxy = (aL - aT)*(v.x*v.y/vlen) + Dm;
        double Dxz = (aL - aT)*(v.x*v.z/vlen) + Dm;

        double Dyx = (aL - aT)*(v.y*v.x/vlen) + Dm;
        double Dyy = aT*vlen + (aL - aT)*(v2y/vlen) + Dm;
        double Dyz = (aL - aT)*(v.y*v.z/vlen) + Dm;

        double Dzx = (aL - aT)*(v.z*v.x/vlen) + Dm;
        double Dzy = (aL - aT)*(v.z*v.y/vlen) + Dm;
        double Dzz = aT*vlen + (aL - aT)*(v2z/vlen) + Dm;
*/

        auto Z1 = static_cast<double>(RG->randomNormal());
        auto Z2 = static_cast<double>(RG->randomNormal());
        auto Z3 = static_cast<double>(RG->randomNormal());
        double xnew = p.x + (v.x + Dxx_p + Dxy_p + Dxz_p) * Dt +
                      (Z1*sqaL * v.x/vlen - Z2 * sqaT * v.y / v2xy - Z3 * sqaT * (v.x * v.z) / (vlen * v2xy)) * sqDTime;
                //std::sqrt(2*Dxx*Dt)*Z1 + std::sqrt(2 * Dxy * Dt) * Z2 + std::sqrt(2 * Dxz * Dt) * Z3;
        double ynew = p.y + (v.y + Dyx_p + Dyy_p + Dyz_p) * Dt +
                      (Z1*sqaL*v.y/vlen + Z2 * sqaT * v.x / v2xy - Z3 * sqaT * (v.y * v.z) / (vlen * v2xy)) * sqDTime;
                      //std::sqrt(2 * Dyx * Dt) * Z1 + std::sqrt(2 * Dyy * Dt) * Z2 + std::sqrt(2 * Dyz * Dt) * Z3;
        double znew = p.z + (v.z + Dzx_p + Dzy_p + Dzz_p) * Dt +
                      (Z1*sqaL*v.z/vlen + Z3 * sqaT * v2xy / vlen) * sqDTime;
                      //std::sqrt(2 * Dzx * Dt) * Z1 + std::sqrt(2 * Dzy * Dt) * Z2 + std::sqrt(2 * Dzz * Dt) * Z3;

        vel.x = xnew - p.x;
        vel.y = ynew - p.y;
        vel.z = znew - p.z;
        pvlu.actualStep = vel.len();
        vel = vel.normalize();
    }
}

#endif //ICHNOS_ICHNOS_RWPT_H
