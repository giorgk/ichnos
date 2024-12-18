#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include <boost/geometry/algorithms/assign.hpp>
#include <boost/math/constants/constants.hpp>

#include "ichnos_structures.h"

#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

namespace ICHNOS {
	namespace DBG {
		/**
		 * Prints the particle positions and velocities are VEX code
		 * \param P
		 * \param prinAttr
		 */
		void displayParticleasVex(const Particle& P, bool prinAttr) {
			vec3 vn = P.getV().normalize();
			std::cout << std::setprecision(5) << std::fixed << "p = addpoint(0,{" << P.getP().x << "," << P.getP().z << "," << P.getP().y << "});";
			if (prinAttr)
				std::cout << std::setprecision(5) << std::fixed << " setpointattrib(0,'N',p,{" << vn.x << "," << vn.z << "," << vn.y << "},'set');";
			std::cout << std::endl;
		}

        /**
         * Prints the particle as VEX code
         * @param p
         */
		void displayVectorasVex(vec3 p) {
			std::cout << std::setprecision(3) << std::fixed << "p = addpoint(0,{" << p.x << "," << p.z << "," << p.y << "});" << std::endl;
		}

        /**
         * Prints the particel and the velocity as VEX code
         * @param p
         * @param v
         */
		void displayPVasVex(vec3 p, vec3 v) {
			vec3 vn = v.normalize();
			std::cout << std::setprecision(3) << std::fixed << "p = addpoint(0,{" << p.x << "," << p.z << "," << p.y << "});";
			std::cout << std::setprecision(3) << std::fixed << " setpointattrib(0,'N',p,{" << vn.x << "," << vn.z << "," << vn.y << "},'set');";
			std::cout << std::endl;
		}

        /**
         * Prints the particle as Matlab code
         * @param p
         */
		void displayVectorMatlab(vec3 p){
		    std::cout << std::setprecision(5) << std::fixed << "[" << p.x << "," << p.y << "," << p.z << "]" << std::endl;
		}
	}

    /**
     * returns the sing of the input number
     * @param num
     * @return
     */
	double sgnFace(int num){
        if (num > 0) return 1.0;
        if (num < 0) return -1.0;
        return 0.0;
	}

    /**
     * Creates a string filled with zeros before the number. Eg for i = 34 and n = 4 this will return 0034
     * @param i Is the number
     * @param n is the total size of the characters.
     * @return
     */
	std::string num2Padstr(int i, int n) {
		std::stringstream ss;
		ss << std::setw(n) << std::setfill('0') << i;
		return ss.str();
	}

    /**
     * Distributes n numbers equally spaced between min and max
     * @param min
     * @param max
     * @param n
     * @param v the distributed numbers
     */
	void linspace(double min, double max, int n, std::vector<double>& v){
		v.clear();
		int iterator = 0;
		for (int i = 0; i <= n-2; ++i){
			double temp = min + i*(max-min)/(floor(static_cast<double>(n)) - 1);
			v.insert(v.begin() + iterator, temp);
        	iterator += 1;
		}
		v.insert(v.begin() + iterator, max);
	}


    /**
     * Checks ith the file is empty
     *
     * See more https://stackoverflow.com/questions/2390912/checking-for-an-empty-file-in-c
     * @param pFile
     * @return true of false
     */
	bool is_file_empty(std::ifstream& pFile){
		return pFile.peek() == std::ifstream::traits_type::eof();
	}

    /**
     * Distributes particles around the well screen layer by layer
     * @param eid The well id
     * @param x x coordinate of the well
     * @param y y coordinate of the well
     * @param top top of the well screen
     * @param bot bottom of the well screen
     * @param S Output vector that holds the particle positions as streamline
     * @param wopt Well options
     * @param releaseTime The release time for the particles
     */
	void distributeParticlesAroundWellLayered(int eid, double x, double y, double top, double bot, std::vector<Streamline>& S, WellOptions wopt, double releaseTime){
		double pi = boost::math::constants::pi<double>();
		std::vector<double> zval;
		linspace(bot, top, wopt.Nlayer, zval);
		int Nppl = wopt.Nparticles/wopt.Nlayer;
		double rads = (2.0*pi)/static_cast<double>(Nppl);
		std::vector<double> rads1; 
		linspace(0.0, 2.0*pi, wopt.Nlayer, rads1);
		std::vector<std::vector<double> > radpos;

		for (unsigned int i = 0; i < static_cast<unsigned int>(wopt.Nlayer); ++i){
			std::vector<double> tmp;
			linspace(0 + rads/2.0 + rads1[i],
								2.0*pi - rads/2.0 + rads1[i], Nppl, tmp);
			radpos.push_back(tmp);
		}
		
		int sid = 0;
		for(int i = 0; i < wopt.Nlayer; ++i){
			for( int j = 0; j < Nppl; ++j){
				vec3 temp;
				temp.x = cos(radpos[i][j])*wopt.Radius + x;
				temp.y = sin(radpos[i][j])*wopt.Radius + y;
				temp.z = zval[i];
				S.push_back(Streamline(eid, sid, Particle(temp, releaseTime)));
				//DEBUG::displayVectorasVex(temp);
				sid = sid + 1;
			}
		}
	}

    /**
     * This is no longer used. Is distributing the particles in a spiral way
     * @param eid
     * @param x
     * @param y
     * @param top
     * @param bot
     * @param S
     * @param wopt
     */
	void distributeParticlesAroundWell(int eid, double x, double y, double top, double bot, std::vector<Streamline>& S, WellOptions wopt) {
		const double PI = 4 * atan(1);
		double maxt = 2*PI*static_cast<double>(wopt.Nlayer);
		double dt = 1 / (static_cast<double>(wopt.Nparticles) - 1);
		double t;
		double r = wopt.Radius;
		double h = top - bot;
		vec3 pp;

		int sid = 0;
		for (double i = 0; i <= 1+dt/2.0; i = i + dt){
			t = i * maxt;
			pp.x = x + r * cos(t);
			pp.y = y + r * sin(t);
			pp.z = bot + h*i;
			//Particle(pp).displayAsVEX(false);
			S.push_back(Streamline(eid, sid, Particle(pp)));
			sid = sid + 1;
		}
	}

    /**
     * Calculates the upper and lower corner of the box that encapsulates the particle.
     * @param p is the particle position
     * @param l is the output lower point
     * @param u is the output upper point
     * @param diameter the size of the box in the horizontal direction
     * @param ratio the ration between horizontal and vertical dimension
     * @param search_mult a multiplier that modifies the dimensions of the bounding box
     */
    void calculate_search_box(vec3 &p, vec3 &l, vec3 &u, double diameter, double ratio, double search_mult) {
        double xy_dir = (diameter/2)*search_mult;
        double z_dir = xy_dir/ratio;
        l.x = p.x - xy_dir;
        l.y = p.y - xy_dir;
        l.z = p.z - z_dir;
        u.x = p.x + xy_dir;
        u.y = p.y + xy_dir;
        u.z = p.z + z_dir;
    }
    /*!
     * Calculates the length of the segment that intersects the bbox along the velocity vector
     * @param p The point that sets the location of the segment in 3D
     * @param v The vector that sets the orientation
     * @param l The lower point of the bounding box
     * @param u The upper point of the bounding box
     * @return the length of the segment
     */
    double diameter_along_velocity(vec3& p, vec3& v, vec3& l, vec3& u){
        // DEBUG::displayPVasVex(p,v);
        // DEBUG::displayVectorasVex(l);
        // DEBUG::displayVectorasVex(u);
        // Find another point along the vector v with origin p
        vec3 p1 = p + v * 0.1;

        double t_max = 1000000000;
        double t_min = -1000000000;
        double tmp_min, tmp_max;

        if (std::abs(v.x) < 0.00000000001){
            tmp_min = -1000000000;
            tmp_max =  1000000000;
        }
        else if (v.x > 0.0 ){
            tmp_min = (l.x - p.x)/(p1.x - p.x);
            tmp_max = (u.x - p.x)/(p1.x - p.x);
        }
        else if (v.x < 0.0){
            tmp_min = (u.x - p.x)/(p1.x - p.x);
            tmp_max = (l.x - p.x)/(p1.x - p.x);
        }
        if (tmp_max < t_max)
            t_max = tmp_max;
        if (tmp_min > t_min)
            t_min = tmp_min;

        if (std::abs(v.y) < 0.00000000001){
            tmp_min = -1000000000;
            tmp_max =  1000000000;
        }
        else if (v.y > 0.0 ){
            tmp_min = (l.y - p.y)/(p1.y - p.y);
            tmp_max = (u.y - p.y)/(p1.y - p.y);
        }
        else if (v.y < 0.0){
            tmp_min = (u.y - p.y)/(p1.y - p.y);
            tmp_max = (l.y - p.y)/(p1.y - p.y);
        }
        if (tmp_max < t_max)
            t_max = tmp_max;
        if (tmp_min > t_min)
            t_min = tmp_min;

        if (std::abs(v.z) < 0.00000000001){
            tmp_min = -1000000000;
            tmp_max =  1000000000;
        }
        else if (v.z > 0.0 ){
            tmp_min = (l.z - p.z)/(p1.z - p.z);
            tmp_max = (u.z - p.z)/(p1.z - p.z);
        }
        else if (v.z < 0.0){
            tmp_min = (u.z - p.z)/(p1.z - p.z);
            tmp_max = (l.z - p.z)/(p1.z - p.z);
        }
        if (tmp_max < t_max)
            t_max = tmp_max;
        if (tmp_min > t_min)
            t_min = tmp_min;

        vec3 pmin = p*(1-t_min) + p1 * t_min;
        vec3 pmax = p*(1-t_max) + p1 * t_max;
        return pmin.distance(pmax.x, pmax.y, pmax.z);
	}

    /**
     * Calculates the barycentric coordinates of a point in triangle. The point can be located outside of the triangle.
     * @param bc The barycentric coordinates
     * @param p The point in question
     * @param p1 First corner of triangle
     * @param p2 Second corner of triangle
     * @param p3 Third corner of triangle
     */
	void calculateBaryCoordTriangle2D(vec3& bc,vec3 &p, vec3 &p1, vec3 &p2, vec3 &p3){
        double det = (p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y);
        if (std::abs(det) < 0.0000001){
            bc.x = -9;
            bc.y = -9;
            bc.z = -9;
            return;
        }
        double bc1 = (p2.y - p3.y)*(p.x - p3.x) + (p3.x - p2.x)*(p.y - p3.y);
        double bc2 = (p3.y - p1.y)*(p.x - p3.x) + (p1.x - p3.x)*(p.y - p3.y);
        bc.x = bc1/det;
        bc.y = bc2/det;
        bc.z = 1 - bc.x - bc.y;
    }

    /**
     * Checks if a point is in or out of a triangle
     * @param p The point in question
     * @param p1 First corner of triangle
     * @param p2 Second corner of triangle
     * @param p3 Third corner of triangle
     * @param bc THe function returs the barycentric coordinates
     * @return true of false
     */
    bool isPointInTriangle(vec3 &p, vec3 &p1, vec3 &p2, vec3 &p3, vec3 &bc){
        calculateBaryCoordTriangle2D(bc, p, p1, p2, p3);
        if (bc.x >= 0 && bc.x <= 1 && bc.y >= 0 && bc.y <= 1 && bc.z >= 0 && bc.z <= 1)
            return true;
        else
            return false;
    }

    /**
     * A helper function which prints information on a console during runtime
     * @param count_times
     * @param FrequencyStat
     * @param calc_time
     * @param max_calc_time
     */
    void PrintStat(int &count_times, int &FrequencyStat, double &calc_time, double &max_calc_time) {
        if (count_times > FrequencyStat){
            std::cout << "||			---Velocity Calc time: " << std::fixed << std::setprecision(15) << calc_time/static_cast<double>(count_times) <<  ", (";// std::endl;
            std::cout << max_calc_time << "), ";
            std::cout << std::endl << std::flush;
            count_times = 0;
            calc_time = 0.0;
            max_calc_time = 0.0;
        }
    }

    /**
     * Calculates the shape functions of a quadrilateral
     * @param u normalized coordinate
     * @param v normalized coordinate
     * @param N1 Shape value
     * @param N2 Shape value
     * @param N3 Shape value
     * @param N4 Shape value
     */
    void QuadShapeFunctions(double u, double v, double &N1, double &N2, double &N3, double &N4){
        N1 = 0.25*(1-u)*(1-v);
        N2 = 0.25*(1+u)*(1-v);
        N3 = 0.25*(1+u)*(1+v);
        N4 = 0.25*(1-u)*(1+v);
    }

    /**
     * Calculates the global coordinate of a point from its local coordinates
     * @param u local coordinate
     * @param v local coordinate
     * @param p1 Corner of quadrilateral
     * @param p2 Corner of quadrilateral
     * @param p3 Corner of quadrilateral
     * @param p4 Corner of quadrilateral
     * @param p Global coordinate
     */
    void local2GlobalCoord(double u, double v, vec3 &p1, vec3 &p2, vec3 &p3, vec3 &p4, vec3 &p){
        p.zero();
        double n1, n2, n3, n4;
        QuadShapeFunctions(u,v, n1, n2, n3, n4);

        p.x = n1*p1.x + n2*p2.x + n3*p3.x + n4*p4.x;
        p.y = n1*p1.y + n2*p2.y + n3*p3.y + n4*p4.y;
    }

    /**
     * Calculates the jacobian of a trangle
     * @param p1
     * @param p2
     * @param p3
     * @param a
     * @param b
     * @param c
     * @param d
     */
    void calculateTriangleJacobian(vec3 &p1, vec3 &p2, vec3 &p3,
                                   double &a, double &b, double &c, double &d){
        a = p1.x - p3.x;
        b = p1.y - p3.y;
        c = p2.x - p3.x;
        d = p2.y - p3.y;
    }

    /**
     * Calculates the jacobian of a quadrilateral
     * @param u
     * @param v
     * @param p1
     * @param p2
     * @param p3
     * @param p4
     * @param a
     * @param b
     * @param c
     * @param d
     */
    void calculateQuadJacobian(double u, double v, vec3 &p1, vec3 &p2, vec3 &p3, vec3 &p4,
                               double &a, double &b, double &c, double &d){
        double dn1 = 0.25*(v - 1);
        double dn2 = 0.25*(1 - v);
        double dn3 = 0.25*(v + 1);
        double dn4 = -0.25*(v + 1);
        double dn5 = 0.25*(u - 1);
        double dn6 = -0.25*(u + 1);
        double dn7 = 0.25*(u + 1);
        double dn8 = 0.25*(1 - u);
        a = dn1*p1.x + dn2*p2.x + dn3*p3.x + dn4*p4.x;
        b = dn1*p1.y + dn2*p2.y + dn3*p3.y + dn4*p4.y;
        c = dn5*p1.x + dn6*p2.x + dn7*p3.x + dn8*p4.x;
        d = dn5*p1.y + dn6*p2.y + dn7*p3.y + dn8*p4.y;
    }

    void calculateQuadInverseJacobian(double u, double v, vec3 &p1, vec3 &p2, vec3 &p3, vec3 &p4,
                                      double &inva, double &invb, double &invc, double &invd){
        double a, b, c, d;
        calculateQuadJacobian(u, v, p1, p2, p3, p4, a, b, c, d);
        double detA = 1/(a*d - b*c);
        inva = d*detA;
        invb = -b*detA;
        invc = -c*detA;
        invd = a*detA;
    }

    bool QuadUVApprox(vec3 &p, vec3 &p1, vec3 &p2, vec3 &p3, vec3 &p4, vec3 &uv){
        bool out = false;
        vec3 p2d(p.x, p.y, 0.0);
        double u = 0;
        double v = 0;
        double u_best = u;
        double v_best = v;
        vec3 puv;
        local2GlobalCoord(u, v, p1, p2, p3, p4, puv);
        double err = p2d.distance(puv.x, puv.y, puv.z);
        double err_best = err;
        double err_prev = err;
        int iter = 0;

        double duv = 0.5;
        std::vector<int> du{-1, -1, 1, 1};
        std::vector<int> dv{-1, 1, -1, 1};
        while (true){
            for (unsigned int i = 0; i < du.size(); ++i){
                u = u_best + du[i]*duv;
                v = v_best + dv[i]*duv;
                local2GlobalCoord(u, v, p1, p2, p3, p4, puv);
                err = p2d.distance(puv.x, puv.y, puv.z);
                if (err < err_best){
                    u_best = u;
                    v_best = v;
                    err_best = err;
                }
            }
            u = u_best;
            v = v_best;
            duv = duv/2;
            if (err_best < 0.5){
                break;
            }
            if (iter > 20 && abs(err_prev - err_best) < 0.01 ){
                break;
            }
            iter = iter + 1;
            if (iter > 50){
                break;
            }
            err_prev = err_best;
        }
        if (err_best < 1){
            uv.x = u_best;
            uv.y = v_best;
            out = true;
        }
        return out;
    }

    int getSign(int a){
        if ( a < 0)
            return -1;
        else
            return 1;
    }

    void QuadInverseMapping(vec3 &p, vec3 &p1, vec3 &p2, vec3 &p3, vec3 &p4, vec3 &uv, double tol = 0.000001){
        // See https://doi.org/10.1115/1.4027667
        vec3 p2d(p.x, p.y, 0.0);
        double u = 0;
        double v = 0;
        double a,b,c,d;
        //DBG::displayVectorMatlab(p1);
        //DBG::displayVectorMatlab(p2);
        //DBG::displayVectorMatlab(p3);
        //DBG::displayVectorMatlab(p4);
        //DBG::displayVectorMatlab(p2d);
        vec3 puv;
        local2GlobalCoord(u, v, p1, p2, p3, p4, puv);
        //DBG::displayVectorMatlab(puv);
        double dst = p2d.distance(puv.x, puv.y, puv.z);
        int iter = 0;
        while (dst > tol) {
            calculateQuadInverseJacobian(u,v,p1,p2,p3,p4,a,b,c,d);
            u = u + a*(p.x - puv.x) + b*(p.y - puv.y);
            v = v + c*(p.x - puv.x) + d*(p.y - puv.y);
            local2GlobalCoord(u, v, p1, p2, p3, p4, puv);
            //DBG::displayVectorMatlab(puv);
            dst = p2d.distance(puv.x, puv.y, puv.z);
            iter = iter + 1;
            if (iter > 2000 || dst > 99999999) {
                bool tf1 = QuadUVApprox(p,p1,p2,p3,p4,uv);
                u = uv.x;
                v = uv.y;
                if (tf1){
                    //std::cout << "QuadUVApprox was used for the inverse parametric mapping" << std::endl;
                }
                else{
                    //std::cout << "Failed to find the parametric coordinates for the point:" << std::endl;
                    //std::cout << std::setprecision(5) << std::fixed << p.x << "," << p.y << std::endl;
                    //std::cout << "for the following quad:" << std::endl;
                    //std::cout << std::setprecision(5) << std::fixed << p1.x << "," << p1.y << std::endl;
                    //std::cout << std::setprecision(5) << std::fixed << p2.x << "," << p2.y << std::endl;
                    //std::cout << std::setprecision(5) << std::fixed << p3.x << "," << p3.y << std::endl;
                    //std::cout << std::setprecision(5) << std::fixed << p4.x << "," << p4.y << std::endl;
                    //std::cout << " uv coords were set to the center of element" << std::endl;
                    u = 0;
                    v = 0;
                }
                break;
            }
        }
        uv.x = u;
        uv.y = v;
    }

    void printFailedInverseQuadInfo(vec3 &p, vec3 &p1, vec3 &p2, vec3 &p3, vec3 &p4){
        std::cout << "Failed to find the parametric coordinates for the point:" << std::endl;
        std::cout << std::setprecision(5) << std::fixed << p.x << "," << p.y << std::endl;
        std::cout << "for the following quad:" << std::endl;
        std::cout << std::setprecision(5) << std::fixed << p1.x << "," << p1.y << std::endl;
        std::cout << std::setprecision(5) << std::fixed << p2.x << "," << p2.y << std::endl;
        std::cout << std::setprecision(5) << std::fixed << p3.x << "," << p3.y << std::endl;
        std::cout << std::setprecision(5) << std::fixed << p4.x << "," << p4.y << std::endl;
        std::cout << " uv coords were set to the center of element" << std::endl;
    }

    void QuadInverseMappingV1(vec3 &p, vec3 &p1, vec3 &p2, vec3 &p3, vec3 &p4, vec3 &uv, bool &tf){
        tf = true;
        //DBG::displayVectorMatlab(p);
        //DBG::displayVectorMatlab(p1);
        //DBG::displayVectorMatlab(p2);
        //DBG::displayVectorMatlab(p3);
        //DBG::displayVectorMatlab(p4);
        double u = 0;
        double v = 0;
        double A = (p1.y - p2.y)*(p3.x - p4.x) - (p3.y - p4.y)*(p1.x - p2.x);
        //A  = (Y(1)-Y(2))*(X(3)-X(4))-(Y(3)-Y(4))*(X(1)-X(2))
        double BO = p1.y*p4.x - p2.y*p3.x + p3.y*p2.x - p4.y*p1.x;
        //BO = Y(1)*X(4)-Y(2)*X(3)+Y(3)*X(2)-Y(4)*X(1)
        double BX = -p1.y + p2.y - p3.y + p4.y;
        //BX = -Y(1)+Y(2)-Y(3)+Y(4)
        double BY = p1.x - p2.x + p3.x - p4.x;
        //BY = X(1)-X(2)+X(3)-X(4)
        double B = BO + BX*p.x + BY*p.y;
        //B  = BO+BX*XP+BY*YP
        double CO = -(p1.y + p2.y)*(p3.x + p4.x) + (p3.y + p4.y)*(p1.x + p2.x);
        //CO = -(Y(1)+Y(2))*(X(3)+X(4))+(Y(3)+Y(4))*(X(1)+X(2))
        double CX = p1.y + p2.y - p3.y - p4.y;
        //CX = Y(1)+Y(2)-Y(3)-Y(4)
        double CY = -p1.x - p2.x + p3.x + p4.x;
        //CY = -X(1)-X(2)+X(3)+X(4)
        double C = CO + 2.0*(CX*p.x + CY*p.y);
        //C  = CO+2.0*(CX*XP+CY*YP)
        if (std::abs(A) < 0.00000001){
            if (std::abs(B) < 0.00000001){
                printFailedInverseQuadInfo(p,p1,p2,p3,p4);
                tf = false;
            }
            else{
                u = -C/(2.0*B);
            }
        }
        else{
            double xt = B*B-A*C;
            if (xt < 0.0){
                printFailedInverseQuadInfo(p,p1,p2,p3,p4);
                tf = false;
            }
            else{
                u = (-B+std::sqrt(xt))/A;
            }
        }
        A = (p2.y - p3.y)*(p1.x - p4.x) - (p1.y - p4.y)*(p2.x - p3.x);
        //A  = (Y(2)-Y(3))*(X(1)-X(4))-(Y(1)-Y(4))*(X(2)-X(3))
        BO = p1.y*p2.x - p2.y*p1.x + p3.y*p4.x - p4.y*p3.x;
        //BO = Y(1)*X(2)-Y(2)*X(1)+Y(3)*X(4)-Y(4)*X(3)
        BX = -p1.y + p2.y - p3.y + p4.y;
        //BX = -Y(1)+Y(2)-Y(3)+Y(4)
        BY = p1.x - p2.x + p3.x - p4.x;
        //BY = X(1)-X(2)+X(3)-X(4)
        B = BO + BX*p.x + BY*p.y;
        //B  = BO+BX*XP+BY*YP
        CO = -(p1.y + p4.y)*(p2.x + p3.x) + (p2.y + p3.y)*(p1.x + p4.x);
        //CO = -(Y(1)+Y(4))*(X(2)+X(3))+(Y(2)+Y(3))*(X(1)+X(4))
        CX = p1.y - p2.y - p3.y + p4.y;
        //CX = Y(1)-Y(2)-Y(3)+Y(4)
        CY = -p1.x + p2.x + p3.x - p4.x;
        //CY = -X(1)+X(2)+X(3)-X(4)
        C = CO + 2.0*(CX*p.x + CY*p.y);
        //C  = CO+2.0*(CX*XP+CY*YP)
        if (std::abs(A) < 0.00000001){
            if (std::abs(B) < 0.00000001){
                printFailedInverseQuadInfo(p,p1,p2,p3,p4);
                tf = false;
            }
            else{
                v = -C/(2.0*B);
            }
        }
        else{
            double yt = B*B-A*C;
            if (yt < 0.0){
                printFailedInverseQuadInfo(p,p1,p2,p3,p4);
                tf = false;
            }
            else{
                v = (-B-std::sqrt(yt))/A;
            }
        }
        uv.x = u;
        uv.y = v;
    }


    std::string getExtension(std::string filename){
        // https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
        // Solution 1.3: stepping away from the iterators
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(filename);
        while (std::getline(tokenStream, token, '.')){
            tokens.push_back(token);
        }
        return tokens.back();
    }

	namespace READ {
		void readTopBot(std::string filename, PointSet2& Pset, bool top, bool bot) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::vector< std::pair<cgal_point_2,elev_data> > xy_data;
				elev_data el_data;
				std::string line;
				double x, y;
				while (getline(datafile, line)) {
					if (line.size() > 1) {
						std::istringstream inp(line.c_str());
						inp >> x;
						inp >> y;
						if (top && bot){
							inp >> el_data.top;
							inp >> el_data.bot;
							if (el_data.top < el_data.bot) {
								std::cerr << "Point:(" << x << "," << y << ") has lower top" << std::endl;
							}
						}
						else if (top){
							inp >> el_data.top;
						}
						else if (bot){
							inp >> el_data.bot;
						}
						cgal_point_2 p(x, y);
						xy_data.push_back(std::make_pair(p, el_data));
					}
				}
				Pset.insert(xy_data.begin(), xy_data.end());
				datafile.close();
			}
		}

		bool readGraphFile(std::string fileGraph, std::vector< std::vector<int> >& data){
            std::cout << "\tReading file " + fileGraph << std::endl;

            std::ifstream datafile(fileGraph.c_str());
            if (!datafile.good()) {
                std::cout << "Can't open the file " << fileGraph << std::endl;
                return false;
            } else {
                std::string line;
                int n, id;
                int idx = 0;
                while (getline(datafile, line)){
                    if (line.size() > 1){
                        std::istringstream inp(line.c_str());
                        inp >> n;
                        std::vector<int>tmp(n+1);
                        tmp[0] = idx;
                        for (int i = 0; i < n; ++i) {
                            inp >> id;
                            tmp[i+1] = id-1;
                        }
                        data.push_back(tmp);
                        idx = idx + 1;
                    }
                }
                return true;
            }
            return false;
		}

		bool readXYZfile(std::string fileXYZ, std::vector<cgal_point_3>& pp, std::vector<Pnt_info>& dd) {
            std::cout << "\tReading file " + fileXYZ << std::endl;


            std::string ext = getExtension(fileXYZ);
            if (ext.compare("h5") == 0){
#if _USEHF > 0
                const std::string XYZDRNameSet("XYZDR");
                const std::string PROCNameSet("PROC");
                HighFive::File HDFNfile(fileXYZ, HighFive::File::ReadOnly);
                HighFive::DataSet datasetXYZDR = HDFNfile.getDataSet(XYZDRNameSet);
                HighFive::DataSet datasetPROC = HDFNfile.getDataSet(PROCNameSet);
                std::vector<std::vector<double>> XYZDR;
                std::vector<int> PROC;
                datasetXYZDR.read(XYZDR);
                datasetPROC.read(PROC);
                Pnt_info td;
                for (int i = 0; i < PROC.size(); ++i) {
                    td.proc = PROC[i];
                    td.diameter = XYZDR[3][i];
                    td.ratio = XYZDR[4][i];
                    td.id = i;
                    pp.push_back(cgal_point_3(XYZDR[0][i], XYZDR[1][i], XYZDR[2][i]));
                    dd.push_back(td);
                }
                return true;
#endif
                return false;
            }


            std::ifstream datafile(fileXYZ.c_str());
            if (!datafile.good()) {
                std::cout << "Can't open the file " << fileXYZ << std::endl;
                return false;
            } else {
                std::string line;
                double x, y, z;
                Pnt_info td;
                int cnt = 0;
                while (getline(datafile, line)) {
                    if (line.size() > 1) {
                        std::istringstream inp(line.c_str());
                        inp >> x;
                        inp >> y;
                        inp >> z;
                        inp >> td.proc;
                        inp >> td.diameter;
                        inp >> td.ratio;
                        td.id = cnt;
                        cnt = cnt + 1;

                        pp.push_back(cgal_point_3(x, y, z));
                        dd.push_back(td);
                    }
                }
                datafile.close();
                return true;
            }
        }

        
		
		/*void read2DscalarField(std::string filename, pointCloud<double>& pntCld) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::string line;
				double x, y, z;
				getline(datafile, line);
				{
					std::istringstream inp(line.c_str());
					inp >> pntCld.Radious;
					inp >> pntCld.Power;
					inp >> pntCld.Threshold;
				}
				while (getline(datafile, line)) {
					std::istringstream inp(line.c_str());
					inp >> x;
					inp >> y;
					inp >> z;
					pntCld.pts.push_back(vec3(x, y, 0.0));
					pntCld.vel.push_back(z);
				}
				datafile.close();
			}
		}*/

		void readNPSATVelocity(std::string filename, 
								std::vector<std::pair<cgal_point, NPSAT_data>>& data, 
								double VelocityMultiplier) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::string line;
				int ii;// , iproc;
				double x, y, z; // , vx, vy, vz;
				NPSAT_data npsat_data;
				ii = 0;
				while (getline(datafile, line)) {
					if (line.size() > 1) {
						std::istringstream inp(line.c_str());
						inp >> x;
						inp >> y;
						inp >> z;
						inp >> npsat_data.v.x;
						inp >> npsat_data.v.y;
						inp >> npsat_data.v.z;
						inp >> npsat_data.proc;
						npsat_data.v = npsat_data.v * VelocityMultiplier;
						npsat_data.id = ii;
						ii = ii + 1;
						data.push_back(std::make_pair(cgal_point(x, y, z), npsat_data));
					}
				}
				std::cout << "Size of Data: " << data.size() << std::endl;
				datafile.close();
			}
		}
		
		void readNPSATVelocity(std::string filename, 
							std::vector<cgal_point_3>& pp,
							std::vector<NPSAT_data>& dd,
								double VelocityMultiplier) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::string line;
				int ii, oldproc;// , iproc;
				double x, y, z; // , vx, vy, vz;
				NPSAT_data npsat_data;
				ii = 0;
				while (getline(datafile, line)) {
					if (line.size() > 1) {
						std::istringstream inp(line.c_str());
						inp >> x;
						inp >> y;
						inp >> z;
						inp >> npsat_data.proc;
						inp >> npsat_data.v.x;
						inp >> npsat_data.v.y;
						inp >> npsat_data.v.z;
						inp >> oldproc;
						//npsat_data.proc = 0;
						inp >> npsat_data.diameter;
						inp >> npsat_data.ratio;
						npsat_data.v = npsat_data.v * VelocityMultiplier;
						npsat_data.id = ii;
						ii = ii + 1;
						//dd.push_back(std::make_pair(cgal_point_3(x, y, z), npsat_data));
						pp.push_back(cgal_point_3(x, y, z));
						dd.push_back(npsat_data);
					}
				}
				std::cout << "Size of Data: " << dd.size() << std::endl;
				datafile.close();
			}
		}
		

		/*void readVelocityFieldFile(std::string filename, pointCloud<vec3>& pntCld) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::string comment("#");
				std::string line;
				int ii;
				double x, y, z;
				int readSection = 0;
				while (getline(datafile, line)) {
					if (line.size() > 1)
						if (comment.compare(0, 1, line, 0, 1) == 0)
							continue;
					if (readSection == 0) {
						std::istringstream inp(line.c_str());
						inp >> pntCld.Radious;
						inp >> pntCld.Power;
						inp >> pntCld.Threshold;
						inp >> ii;
						inp >> pntCld.NmaxPnts;
						inp >> pntCld.NminPnts;
						if (ii == 0)
							pntCld.InterpolateOutside = false;
						else
							pntCld.InterpolateOutside = true;
						readSection++;
						continue;
					}

					if (readSection == 1) {
						std::istringstream inp(line.c_str());
						inp >> ii;
						inp >> x;
						inp >> y;
						inp >> z;
						pntCld.proc.push_back(ii);
						pntCld.pts.push_back(vec3(x, y, z));
						inp >> x;
						inp >> y;
						inp >> z;
						pntCld.vel.push_back(vec3(x, y, z));
					}
				}
				datafile.close();

			}
		}*/

		void readProcessorDomain(std::string filename, boostPolygon& polygon, int myRank) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the ProcessorPolys file:[" << filename << "]" << std::endl;
			}
			else {
				std::vector<boostPoint> polygonPoints;
				std::string line;
				getline(datafile, line);
				int Npoly, id, nVerts;
				double x, y;
				{
					std::istringstream inp(line.c_str());
					inp >> Npoly;
				}
				bool myPolyFound = false;
				for (int i = 0; i < Npoly; ++i) {
					getline(datafile, line);
					{
						std::istringstream inp(line.c_str());
						inp >> id;
						inp >> nVerts;
						if (id == myRank){
                            myPolyFound = true;
						}
					}
					for (int j = 0; j < nVerts; ++j) {
						getline(datafile, line);
						if (myPolyFound) {
							std::istringstream inp(line.c_str());
							inp >> x;
							inp >> y;
							polygonPoints.push_back(boostPoint(x, y));
						}
					}
					if (myPolyFound) {
						boost::geometry::assign_points(polygon, polygonPoints);
                        //std::cout << "Area Before" << boost::geometry::area(polygon) << std::endl;
						boost::geometry::correct(polygon);
                        //std::cout << "Area After" << boost::geometry::area(polygon) << std::endl;
                        //bool tf = boost::geometry::within(boostPoint(-30818.557, -183282.562), polygon);
                        //std::cout << myRank << " " << tf << std::endl;
                        //if (tf){
                        //    std::cout << polygonPoints[0].x() << "," << polygonPoints[0].y() << std::endl;
                        //    std::cout << polygonPoints[1].x() << "," << polygonPoints[1].y() << std::endl;
                        //    std::cout << polygonPoints[2].x() << "," << polygonPoints[2].y() << std::endl;
                        //    std::cout << polygonPoints[3].x() << "," << polygonPoints[3].y() << std::endl;
                        //}
						break;
					}
				}
				datafile.close();
			}
		}


		bool readParticleFile(std::string filename, std::vector <Streamline>& S, bool isTransient) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file " << filename << std::endl;
				return false;
			}
			else {
				std::string comment("#");
				std::string line;
				int eid, sid;
				double rt = 0;
				vec3 pp;
				while (getline(datafile, line)) {
					if (line.size() >= 1 && !line.empty())
						if (comment.compare(0, 1, line, 0, 1) == 0)
							continue;
                    if (line.empty())
                        continue;

					std::istringstream inp(line.c_str());
					inp >> eid;
					inp >> sid;
					inp >> pp.x;
					inp >> pp.y;
					inp >> pp.z;
					if (isTransient)
					    inp >> rt;
					//P.push_back(Particle<T>(pp));
					S.push_back(Streamline(eid, sid, Particle(pp, rt)));
				}
				datafile.close();
				return true;
			}
		}

		bool readWellFile(std::string filename, std::vector <Streamline>& S, bool isTransient) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file " << filename << std::endl;
				return false;
			}
			else {
                std::string comment("#");
				WellOptions wopt;
				std::string line;
				int eid;
				double rt = 0;
				double x, y, top, bot;
				while (getline(datafile, line)){
                    if (line.size() >= 1 && !line.empty())
                        if (comment.compare(0, 1, line, 0, 1) == 0)
                            continue;
                    if (line.empty())
                        continue;
					std::istringstream inp(line.c_str());
					inp >> wopt.Nparticles;
					inp >> wopt.Nlayer;
					inp >> wopt.Radius;
					break;
				}
				while (getline(datafile, line)) {
                    if (line.size() >= 1 && !line.empty())
                        if (comment.compare(0, 1, line, 0, 1) == 0)
                            continue;
                    if (line.empty())
                        continue;
					std::istringstream inp(line.c_str());
					inp >> eid;
					inp >> x;
					inp >> y;
					inp >> top;
					inp >> bot;
                    if (isTransient)
                        inp >> rt;
					distributeParticlesAroundWellLayered(eid, x, y, top, bot, S, wopt, rt);
					//distributeParticlesAroundWell(eid, x, y, top, bot, S, wopt);
				}
				datafile.close();
				return true;
			}
		}

		bool readTrajfiles(std::string filename, streamlineMap& Smap) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
				return false;
			}
			else {
				if (is_file_empty(datafile)){
					return true;
				}
				std::string line, er_str;
				int eid, sid, pid, eid_prev, sid_prev, pid_prev;
				gPart particle;
				gStream streamline;
				streamlineMap::iterator eIt;
				std::map<int, gStream >::iterator sIt;
				eid_prev = -9;
				sid_prev = -9;
				while (getline(datafile, line)) {
					std::istringstream inp(line.c_str());
					inp >> pid;
					inp >> eid;
					inp >> sid;
					if (pid < 0) {
						inp >> er_str;
					}
					else {
                        pid_prev = pid;
						inp >> particle.p.x;
						inp >> particle.p.y;
						inp >> particle.p.z;
						inp >> particle.v.x;
						inp >> particle.v.y;
						inp >> particle.v.z;
						inp >> particle.age;
					}

					bool dofind = true;
					if (eid_prev == eid && sid_prev == sid)
						dofind = false;

					if (dofind)
						eIt = Smap.find(eid);

					if (eIt != Smap.end()) {
						if (dofind)
							sIt = eIt->second.find(sid);
						if (sIt != eIt->second.end()) {
							if (pid < 0) {
                                int lastPID = sIt->second.particles.rbegin()->first;
                                if (lastPID == pid_prev)
								    sIt->second.ex = castExitReasons2Enum(er_str);
							}
							else {
								sIt->second.particles.insert(std::pair<int, gPart>(pid, particle));
							}
						}
						else {
							gStream gs;
							if (pid < 0) {
								gs.ex = castExitReasons2Enum(er_str);
							}
							else {
								gs.particles.insert(std::pair<int, gPart>(pid, particle));
							}
							eIt->second.insert(std::pair<int, gStream >(sid, gs));
							sIt = eIt->second.find(sid);
						}
					}
					else {
						gStream gs;
						if (pid < 0) {
							gs.ex = castExitReasons2Enum(er_str);
						}
						else {
							gs.particles.insert(std::pair<int, gPart>(pid, particle));
						}
						std::map<int, gStream > tmpStream;
						tmpStream.insert(std::pair<int, gStream>(sid, gs));
						Smap.insert(std::pair<int, std::map<int, gStream > >(eid, tmpStream));
						// make sure that the iterators are pointing to the new data that were just inserted
						// most likely the following lines will belong to the same streamline
						eIt = Smap.find(eid);
						sIt = eIt->second.find(sid);
					}
					sid_prev = sid;
					eid_prev = eid;
				}
				return true;
			}
		}

		template<typename T>
		bool read2Darray(std::string filename, int ncols, std::vector<std::vector<T>> &data){
            std::ifstream datafile(filename.c_str());
            if (!datafile.good()) {
                std::cout << "Can't open the file " << filename << std::endl;
                return false;
            }
            else{
                T d;
                data.clear();
                std::vector<T> vecd(ncols);
                std::string line;
                while (getline(datafile, line)){
                    if (line.size() > 1){
                        std::istringstream inp(line.c_str());
                        for (int i = 0; i < ncols; ++i) {
                            inp >> d;
                            vecd[i] = d;
                        }
                        data.push_back(vecd);
                    }
                }
            }
            return true;
		}

		bool readTimeStepFile(std::string filename, std::vector<double>& TS){
            std::ifstream datafile(filename.c_str());
            if (!datafile.good()){
                std::cout << "Can't open the file " << filename << std::endl;
                return false;
            }
            else{
                std::string line;
                double v;
                while (getline(datafile, line)){
                    if (line.size() > 0){
                        std::istringstream inp(line.c_str());
                        inp >> v;
                        TS.push_back(v);
                    }
                }
                datafile.close();
            }
            return true;
		}

		bool read1DVelocity(std::string filename, /*int nPoints,*/ int nSteps, double multiplier, coordDim dim, VelTR& VEL){
            std::vector<std::vector<double>> data;
            bool tf = read2Darray<double>(filename, nSteps, data);
            if (tf){
                int nPoints = static_cast<int>(data.size());
                for (int i = 0; i < nPoints; ++i){
                    for (int j = 0; j < nSteps; ++j){
                        VEL.setVELvalue(data[i][j]*multiplier, i, j, dim);
                    }
                }
                return true;
            }
            return false;
            /*
            std::cout << "\tReading file " + filename << std::endl;
            std::ifstream datafile(filename.c_str());
            if (!datafile.good()){
                std::cout << "Can't open the file " << filename << std::endl;
                return false;
            }
            else{
                std::string line;
                double v;
                for (int i = 0; i < nPoints; i++){
                    getline(datafile, line);
                    std::istringstream inp(line.c_str());
                    for (int j = 0; j < nSteps; j++){
                        inp >> v;
                        VEL.setVELvalue(v*multiplier, i, j, dim);
                    }
                }
            }
            return true;
            */
		}

		bool H5SteadyState3DVelocity(std::string filename, std::vector<std::vector<double>>& VXYZ){
            std::cout << "\tReading file " + filename << std::endl;
#if _USEHF > 0
            const std::string VXYZNameSet("VXYZ");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetVXYZ = HDFNfile.getDataSet(VXYZNameSet);
            //std::vector<std::vector<double>> VXYZ;
            datasetVXYZ.read(VXYZ);
            //for (int i = 0; i < VXYZ[0].size(); i++){
            //    VEL.setVELvalue(VXYZ[0][i]*multiplier, i, 0, coordDim::vx);
            //    VEL.setVELvalue(VXYZ[1][i]*multiplier, i, 0, coordDim::vy);
            //    VEL.setVELvalue(VXYZ[2][i]*multiplier, i, 0, coordDim::vz);
            //}
            return true;
#endif
            return false;
		}

		bool ASCIISteadState3DVelocity(std::string filename, int nPoints, int nInfoSkip, double multiplier, VelTR& VEL){
            std::cout << "\tReading file " + filename << std::endl;
            std::ifstream datafile(filename.c_str());
            if (!datafile.good()){
                std::cout << "Can't open the file " << filename << std::endl;
                return false;
            }
            else{
                std::string line;
                double tmp, vx, vy, vz;
                for (int i = 0; i < nPoints; ++i){
                    getline(datafile, line);
                    std::istringstream inp(line.c_str());
                    for (int j = 0; j < nInfoSkip; ++j){
                        inp >> tmp;
                    }
                    inp >> vx;
                    inp >> vy;
                    inp >> vz;
                    VEL.setVELvalue(vx*multiplier, i, 0, coordDim::vx);
                    VEL.setVELvalue(vy*multiplier, i, 0, coordDim::vy);
                    VEL.setVELvalue(vz*multiplier, i, 0, coordDim::vz);
                }
            }
            datafile.close();
            return true;
		}

		bool H5Transient3DVelocity(std::string filename, /*int nPoints,*/ int nSteps, double multiplier, VelTR& VEL){
            std::cout << "\tReading file " + filename << std::endl;
#if _USEHF > 0
            //VEL.init(nPoints,nSteps);
            const std::string VXNameSet("VX");
            const std::string VYNameSet("VY");
            const std::string VZNameSet("VZ");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetVX = HDFNfile.getDataSet(VXNameSet);
            HighFive::DataSet datasetVY = HDFNfile.getDataSet(VYNameSet);
            HighFive::DataSet datasetVZ = HDFNfile.getDataSet(VZNameSet);
            std::vector<std::vector<double>> VX;
            std::vector<std::vector<double>> VY;
            std::vector<std::vector<double>> VZ;
            datasetVX.read(VX);
            datasetVY.read(VY);
            datasetVZ.read(VZ);


            VEL.init(VX[0].size(), nSteps, 3);

            if (VX.size() != nSteps){
                std::cout << "The number of time steps (" << VX.size() << ") in the file " << filename
                << "\n does not match the number of steps in the Time step file " << nSteps << std::endl;
                return false;
            }
            for (int i = 0; i < VX[0].size(); ++i){
                for (int j = 0; j < VX.size(); ++j) {
                    VEL.setVELvalue(VX[j][i]*multiplier, i, j, ICHNOS::coordDim::vx);
                    VEL.setVELvalue(VY[j][i]*multiplier, i, j, ICHNOS::coordDim::vy);
                    VEL.setVELvalue(VZ[j][i]*multiplier, i, j, ICHNOS::coordDim::vz);
                }
            }
            return true;
#endif
            return false;
		}

	}

	

	namespace WRITE {
		void PrintParticle2Log(std::ofstream& log_file, Streamline& S, int i, PrintOptions printOpt) {
			Particle pp = S.getParticle(i);
			log_file << pp.getPid() << " "
				<< S.getEid() << " "
				<< S.getSid() << " "
				<< std::setprecision(printOpt.Pprec) << std::fixed
				<< pp.getP().x << " " << pp.getP().y << " " << pp.getP().z << " "
				<< std::setprecision(printOpt.Vprec) << std::scientific
				<< pp.getV().x << " " << pp.getV().y << " " << pp.getV().z << " "
                //<< pp.getV().len() << " "
                << std::setprecision(printOpt.Tprec) << std::scientific
				<< pp.getTime() << std::endl;
		}
		void PrintExitReason(std::ofstream& log_file, Streamline& S, ExitReason er) {
			log_file << -9 << " "
				<< S.getEid() << " "
				<< S.getSid() << " " 
				<< castExitReasons2String(er) << std::endl;
		}

		void printStreamslineMap(std::string filename, streamlineMap& SM) {
			std::ofstream out_file;
			out_file.open(filename.c_str());
			streamlineMap::iterator eit = SM.begin();
			std::map<int, gStream >::iterator sit;
			std::map<int, gPart>::iterator pit;
			
			for (; eit != SM.end(); ++eit) {
				sit = eit->second.begin();
				for (; sit != eit->second.end(); ++sit) {
					pit = sit->second.particles.begin();
					//double age = 0;
					//int i = 0;
					//vec3 p_prev, p_curr;
					//vec3 v_prev, v_curr, v_m;
					for (; pit != sit->second.particles.end(); ++pit) {
						/* // This is used to calculate the age
						// However to save some space in the output files
						// I wont print it
						if (i == 0) {
							p_prev = pit->second.p;
							v_prev = pit->second.v;
						}
						else {
							p_curr = pit->second.p;
							v_curr = pit->second.v;
							v_m = (v_prev + v_curr) * 0.5;
							double dst = (p_curr - p_prev).len();
							age += dst / v_m.len();
							p_prev = p_curr;
							v_prev = v_curr;
						}
						*/
						out_file << eit->first << " "
							<< sit->first << " "
							<< std::setprecision(2) << std::fixed
							<< pit->second.p.x << " "
							<< pit->second.p.y << " "
							<< pit->second.p.z << " "
							<< std::setprecision(5) << std::fixed
							<< pit->second.v.len() << " "
							<< pit->second.age
							// << pit->second.v.x << " "
							// << pit->second.v.y << " "
							// << pit->second.v.z << " "
							/* << age */ << std::endl;
						//i++;
					}
					out_file << "-9 " 
							 << eit->first << " " 
							 << sit->first << " "
							 << castExitReasons2String(sit->second.ex)
							 << std::endl;
				}
			}
			out_file.close();
		}

		/**
		 *
		 * @param S
		 * @param filename The name of the file without the extension
		 */
		void writeStreamlines(std::vector<Streamline>& S,
                              std::string filename,
                              PrintOptions printOpt,
                              bool append = false){
            //std::cout << "S size = " << S.size() << std::endl;
#if _USEHF > 0
            if (printOpt.printH5){ //print with Highfive
                const std::string FILE_NAME(filename + ".h5");
                const std::string DATASET_NAME("PVA");
                const std::string DATASET_IDS("ESID");
                const std::string DATASET_EXIT("EXIT");
                std::vector<std::vector<float>> alldata;
                std::vector<std::vector<int>> ids;
                std::vector<std::string> ExitReason;
                //boost::multi_array<double, 2> my_array(boost::extents[size_x][size_y]);
                for (unsigned int i = 0; i < S.size(); ++i){
                    if (!S[i].printIt())
                        continue;
                    std::vector<float> tmp_data(10,0);
                    std::vector<int> tmp_ids(2,0);
                    for (unsigned int j = 0; j < S[i].size(); ++j){
                        tmp_data[0] = static_cast<float>(S[i].getEid());
                        tmp_data[1] = static_cast<float>(S[i].getSid());
                        Particle pp = S[i].getParticle(j);
                        tmp_data[2] = static_cast<float>(pp.getPid());
                        tmp_data[3] = static_cast<float>(pp.getP().x);
                        tmp_data[4] = static_cast<float>(pp.getP().y);
                        tmp_data[5] = static_cast<float>(pp.getP().z);
                        tmp_data[6] = static_cast<float>(pp.getV().x);
                        tmp_data[7] = static_cast<float>(pp.getV().y);
                        tmp_data[8] = static_cast<float>(pp.getV().z);
                        tmp_data[9] = static_cast<float>(pp.getTime());
                        alldata.push_back(tmp_data);
                    }
                    tmp_ids[0] = S[i].getEid();
                    tmp_ids[1] = S[i].getSid();
                    ids.push_back(tmp_ids);
                    ExitReason.push_back(castExitReasons2String(S[i].getExitReason()));
                }
                HighFive::File file(FILE_NAME, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
                HighFive::DataSet dataset = file.createDataSet<float>(DATASET_NAME, HighFive::DataSpace::From(alldata));
                HighFive::DataSet datasetIDS = file.createDataSet<int>(DATASET_IDS, HighFive::DataSpace::From(ids));
                HighFive::DataSet datasetEXIT = file.createDataSet<std::string>(DATASET_EXIT, HighFive::DataSpace::From(ExitReason));
                dataset.write(alldata);
                datasetIDS.write(ids);
                datasetEXIT.write(ExitReason);
            }
#endif
            if (printOpt.printASCII){
                const std::string log_file_name = (filename + ".traj");
                //std::cout << "Printing Output to " << log_file_name << std::endl;
                std::ofstream log_file;
                if (append){
                    log_file.open(log_file_name.c_str(),std::ios_base::app);
                }
                else{
                    log_file.open(log_file_name.c_str());
                }

                for (unsigned int i = 0; i < S.size(); ++i){
                    if (!S[i].printIt())
                        continue;
                    for (unsigned int j = 0; j < S[i].size()-1; ++j) {
                        WRITE::PrintParticle2Log(log_file, S[i], j, printOpt);
                    }
                    WRITE::PrintParticle2Log(log_file, S[i], S[i].size() - 1, printOpt);
                    WRITE::PrintExitReason(log_file, S[i], S[i].getExitReason());
                }
                log_file.close();
            }
		}
	}

	void gather_particles(ICHNOS::options& opt) {
		for (int ireal = 0; ireal < opt.Popt.Nrealizations; ++ireal) {
			streamlineMap Smap;
			for (int i = 0; i < opt.niter; ++i) {
				for (int j = 0; j < opt.nproc; ++j) {
					std::string filename = opt.Popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) +
						"_iter_" + num2Padstr(i, 4) + "_proc_" + num2Padstr(j, 4) + ".traj";

					std::cout << "Reading ... " << filename << std::endl;
					bool tf = READ::readTrajfiles(filename, Smap);
					if (!tf) {
						std::cout << "Error while reading trajectory files" << std::endl;
						return;
					}
				}
				if (!opt.Popt.gatherOneFile) {
					std::string filename = opt.Popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) + "_gather_iter_" + num2Padstr(i, 4) + ".traj";
					ICHNOS::WRITE::printStreamslineMap(filename, Smap);
					Smap.clear();
				}
			}
			if (opt.Popt.gatherOneFile) {
				std::string filename = opt.Popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) + "_gather_ALL.traj";
				ICHNOS::WRITE::printStreamslineMap(filename, Smap);
			}
		}
	}

	//double interpolateScatter2D(cgal_Delaunay_2& T, coord_map& values, double x, double y) {
	//	cgal_point_2 p(x, y);
	//	std::vector<std::pair<cgal_point_2, coord_type>> coords;
	//	coord_type norm = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords)).second;
	//	coord_type res = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, value_access(values));
	//	return static_cast<double>(res);
	//}
	
	/*double interpolateScalarTree(std::unique_ptr<nano_kd_tree_scalar>& tree, vec3 p) {
		double thres = tree->dataset.Threshold;
		double radius = tree->dataset.Radious;
		double power = tree->dataset.Power;
		double query_pt[3] = { p.x, p.y, 0.0 };
		std::vector<std::pair<size_t, double> >   ret_matches;
		nanoflann::SearchParams params;
		size_t N = tree->radiusSearch(&query_pt[0], radius*radius, ret_matches, params);

		if (N == 0) {
			std::cerr << "There are no points around (" << p.x << "," << p.y << ") whithin " << radius << "distance" << std::endl;
			std::cerr << "consider increasing the radius" << std::endl;
			return -999999;
		}
		else {
			double sumW = 0;
			double sumWVal = 0;
			double w;
			for (size_t i = 0; i < N; i++) {
				if (ret_matches[i].second < thres) {
					return tree->dataset.kdtree_get_pt(ret_matches[i].first, 3);
				}
				else {
					w = 1 / std::pow(std::sqrt(ret_matches[i].second), power);
					sumW += w;
					sumWVal += tree->dataset.kdtree_get_vel_vec(ret_matches[i].first) * w;
				}
			}
			return(sumWVal / sumW);
		}
	}*/

	//void interpolateIntTree(std::map<int, double>& id_dst, std::map<int, double>& proc_map,
	//	std::unique_ptr <nano_kd_tree_int>& tree, vec3 p) {
	//	double power = tree->dataset.Power;
	//	double threshold = tree->dataset.Threshold;
	//	double radius = tree->dataset.Radious;
	//	id_dst.clear();
	//	proc_map.clear();
	//	double query_pt[3] = { p.x, p.y, 0.0 };
	//	// This is a vector of pairs (id distance).  
	//	std::vector<std::pair<size_t, double> >   ret_matches;
	//	nanoflann::SearchParams params;
	//	size_t N = tree->radiusSearch(&query_pt[0], radius * radius, ret_matches, params);
	//	if (N == 0) {
	//		std::cerr << "There are no points around (" << p.x << "," << p.y << ") whithin " << radius << "distance" << std::endl;
	//		std::cerr << "consider increasing the radius" << std::endl;
	//	}
	//	else {
	//		double sumW = 0.0;
	//		double w;
	//		int iproc;
	//		std::map<int, double>::iterator it;
	//		for (size_t i = 0; i < N; i++) {
	//			iproc = tree->dataset.kdtree_get_proc(ret_matches[i].first);
	//			if (ret_matches[i].second < threshold) {
	//				id_dst.clear();
	//				id_dst.insert(std::pair<int, double>(static_cast<int>(ret_matches[i].first), 1.0));
	//				break;
	//			}
	//			else {
	//				w = 1 / std::pow(std::sqrt(ret_matches[i].second), power);
	//				it = proc_map.find(iproc);
	//				if (it == proc_map.end()) {
	//					proc_map.insert(std::pair<int, double>(iproc, w));
	//				}
	//				else {
	//					it->second += w;
	//				}
	//				sumW += w;
	//				id_dst.insert(std::pair<int, double>(static_cast<int>(ret_matches[i].first), w));
	//			}
	//		}
	//		it = id_dst.begin();
	//		for (it; it != id_dst.end(); ++it) {
	//			it->second = it->second / sumW;
	//		}
	//		it = proc_map.begin();
	//		for (it; it != proc_map.end(); ++it) {
	//			it->second = it->second / sumW;
	//		}
	//	}
	//}

	//void interpolateVectorTree(vec3& vel, std::map<int, double>& proc_map,
	//	std::unique_ptr <nano_kd_tree_vector>& tree, vec3 p, float scale = 1) {
	//	double power = tree->dataset.Power;
	//	double threshold = tree->dataset.Threshold;
	//	double radius = tree->dataset.Radious;
	//	size_t Nmax = static_cast<size_t>(tree->dataset.NmaxPnts);
	//	size_t Nmin = static_cast<size_t>(tree->dataset.NminPnts);
	//	proc_map.clear();
	//	double query_pt[3] = { p.x, p.y, p.z };

	//	std::vector<std::pair<size_t, double> >   ret_matches;
	//	nanoflann::SearchParams params;
	//	size_t N = 0;
	//	int count_attempts = 0;
	//	while (N < Nmin) {
	//		N = tree->radiusSearch(&query_pt[0], radius * radius, ret_matches, params);
	//		radius = radius * 2;
	//		count_attempts++;
	//		if (count_attempts > 20) {
	//			std::cout << "There are no points around (" << p.x << "," << p.y << "," << p.z << ") whithin " << radius << "  units radius" << std::endl;
	//			vel = vec3();
	//			return;
	//		}
	//	}
	//	N = std::min(N, Nmax);

	//	// find the largest horizontal and vertical distance between the point in question and the
	//	// closest points around this
	//	std::vector<vec3> closest_pnts(N);
	//	double xymax = 0;
	//	double zmax = 0;
	//	double tmpxy, tmpz;
	//	for (size_t i = 0; i < N; i++) {
	//		//std::cout << tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).x << ", "
	//		//	<< tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).y << ", "
	//		//	<< tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).z << std::endl;
	//		closest_pnts[i] = tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first) - p;
	//		tmpxy = std::sqrt(closest_pnts[i].x * closest_pnts[i].x + closest_pnts[i].y * closest_pnts[i].y);
	//		tmpz = std::sqrt(closest_pnts[i].z * closest_pnts[i].z);
	//		if (tmpxy > xymax)
	//			xymax = tmpxy;
	//		if (tmpz > zmax)
	//			zmax = tmpz;
	//	}

	//	double ratio = xymax / zmax;

	//	// interpolate using the scaled interpolate
	//	double w, d;
	//	double sumW = 0;
	//	vec3 sumWVal;

	//	std::map<int, double>::iterator it;
	//	int iproc;

	//	for (size_t i = 0; i < N; i++) {
	//		closest_pnts[i].z = closest_pnts[i].z * ratio * scale + closest_pnts[i].z * (1.0 - scale);
	//		d = std::sqrt(closest_pnts[i].x * closest_pnts[i].x +
	//			closest_pnts[i].y * closest_pnts[i].y +
	//			(closest_pnts[i].z * closest_pnts[i].z));
	//		iproc = tree->dataset.kdtree_get_proc(ret_matches[i].first);
	//		if (d < threshold) {
	//			vel = tree->dataset.kdtree_get_vel_vec(ret_matches[i].first);
	//		}
	//		else {
	//			w = 1 / std::pow(d, power);
	//			it = proc_map.find(iproc);
	//			if (it == proc_map.end()) {
	//				proc_map.insert(std::pair<int, double>(iproc, w));
	//			}
	//			else {
	//				it->second += w;
	//			}
	//			sumW += w;
	//			sumWVal = sumWVal + tree->dataset.kdtree_get_vel_vec(ret_matches[i].first) * w;
	//		}
	//	}

	//	it = proc_map.begin();
	//	for (it; it != proc_map.end(); ++it) {
	//		it->second = it->second / sumW;
	//	}

	//	vel = sumWVal * (1 / sumW);
	//	
	//}

	bool is_input_scalar(std::string input) {
		// try to convert the input to scalar
		bool outcome;
		try {
			double value = std::stod(input);
			value = value + 1; // something to surpress the warning
			outcome = true;
		}
		catch (...) {
			outcome = false;
		}
		return outcome;
	}

	namespace MPI {
		namespace DBG {
			template <typename T>
			void DataPerProc(std::vector<T>& v, int myrank, std::string vname) {
				for (unsigned int i = 0; i < v.size(); ++i)
					std::cout << "Proc: " << myrank << ": " << vname << "[" << i << "] = " << v[i] << std::endl;
			}

			template <typename T>
			void PrintVectorData(std::vector<std::vector<T> >& v, int myrank, std::string vname) {
				for (int i = 0; i < v.size(); ++i){
					for (int j = 0; j < v[i].size(); ++j) {
						std::cout << "Proc: " << myrank << ": " << vname << "[" << i << "][" << j << "] = " << v[i][j] << std::endl;
					}
				}
			}

			template<typename T>
			void PrintVectorSize(std::vector<std::vector<T> >& v, int myrank, std::string vname) {
				if (myrank < static_cast<int>(v.size()))
					std::cout << "I'm proc " << myrank << " and I have " << v[myrank].size() << " elements of " << vname << std::endl;
				else
					std::cout << "Vector " << vname << " doesnt have space for rank " << myrank << std::endl;
			}
		}
		

		

		void Send_receive_size(int N, int n_proc, std::vector<int>& output, boost::mpi::communicator& world) {

			output.clear();
			output.resize(n_proc);
			std::vector<int> temp(n_proc, 1);
			std::vector<int> displs(n_proc);
			for (int i = 1; i < n_proc; ++i)
				displs[i] = displs[i - 1] + 1;

			MPI_Allgatherv(&N, // This is what this processor will send to every other
				1, //This is the size of the message from this processor
				MPI_INT, // The data type will be sent
				&output[0], // This is where the data will be send on each processor
				&temp[0], // an array with the number of points to be sent/receive
				&displs[0],
				MPI_INT, // The data type will be received
				world);

			//for (int i = 0; i < n_proc; ++i)
			//	std::cout << world.rank() << ":" << output[i] << std::endl;
		}

		template <typename T>
		void Send_receive_data(std::vector<std::vector<T> >& data,
			std::vector <int> N_data_per_proc,
			unsigned int my_rank,
			boost::mpi::communicator& world,
			MPI_Datatype MPI_TYPE) {

			// data is a vector of vectors of type T with size equal to n_proc.
			// This function transfer to all processors the content of data[my_rank]
			// if there are any data in data[i], where i=[1,n_proc; i!=myrank] this will be deleted
			// The size of data[my_rank].size() = N_data_per_proc[my_rank]. This is user responsibility

			int N = data[my_rank].size();
			if (N == 0)
				data[my_rank].push_back(0);

			unsigned int n_proc = data.size();
			std::vector<int> displs(n_proc);
			displs[0] = 0;
			for (unsigned int i = 1; i < n_proc; ++i)
				displs[i] = displs[i - 1] + N_data_per_proc[i - 1];
			int totdata = displs[n_proc - 1] + N_data_per_proc[n_proc - 1];
			
			std::vector<T> temp_receive(totdata);

			
			MPI_Allgatherv(&data[my_rank][0], // This is what this processor will send to every other
				N, //This is the size of the message from this processor
				MPI_TYPE, // The data type will be sent
				&temp_receive[0], // This is where the data will be send on each processor
				&N_data_per_proc[0], // an array with the number of points to be sent/receive
				&displs[0],
				MPI_TYPE, world);

			//for (int i = 0; i < temp_receive.size(); ++i)
			//	std::cout << world.rank() << ": temp_receive[" << i << "] = " << temp_receive[i] << std::endl;

			
			// Now put the data in the data vector
			for (unsigned int i = 0; i < n_proc; ++i) {
				data[i].clear();
				data[i].resize(N_data_per_proc[i]);
				for (int j = 0; j < N_data_per_proc[i]; ++j)
					data[i][j] = temp_receive[displs[i] + j];
			}
		}

		template <typename T>
		void broadcast_vector (std::vector<T> & v, int proc, boost::mpi::communicator& world){
			T value;
			int VectorSize;

			if (world.rank() == proc) 
				VectorSize = static_cast<int>(v.size());

			boost::mpi::broadcast(world, VectorSize, proc);
			if (world.rank() != proc)
				v.clear();

			for (int i = 0; i < VectorSize; ++i) {
				if (world.rank() == proc) {
					value = v[i];
				}
				boost::mpi::broadcast(world, value, proc);
				if (world.rank() != proc) {
					v.push_back(value);
				}
			}
		}

		template <typename T>
		void sumScalar(T& scalar, int n_proc, boost::mpi::communicator& world, MPI_Datatype MPI_TYPE) {
			std::vector<T> Allscalar(n_proc);
			std::vector<int> temp(n_proc, 1);
			std::vector<int> displs(n_proc);
			for (int i = 1; i < n_proc; ++i)
				displs[i] = displs[i - 1] + 1;

			MPI_Allgatherv(&scalar, // This is what this processor will send to every other
				1, //This is the size of the message from this processor
				MPI_TYPE, // The data type will be sent
				&Allscalar[0], // This is where the data will be send on each processor
				&temp[0], // an array with the number of points to be sent/receive
				&displs[0],
				MPI_TYPE, // The data type will be received
				world);
			scalar = 0;
			for (int i = 0; i < n_proc; ++i) {
				scalar += Allscalar[i];
			}
		}

		template <typename T>
		void maxScalar(T& scalar, int n_proc, boost::mpi::communicator& world, MPI_Datatype MPI_TYPE) {
			std::vector<T> Allscalar(n_proc);
			std::vector<int> temp(n_proc, 1);
			std::vector<int> displs(n_proc);
			for (unsigned int i = 1; i < n_proc; ++i)
				displs[i] = displs[i - 1] + 1;

			MPI_Allgatherv(&scalar, // This is what this processor will send to every other
				1, //This is the size of the message from this processor
				MPI_TYPE, // The data type will be sent
				&Allscalar[0], // This is where the data will be send on each processor
				&temp[0], // an array with the number of points to be sent/receive
				&displs[0],
				MPI_TYPE, // The data type will be received
				world);
			scalar = -9999999999;
			for (int i = 0; i < n_proc; ++i) {
				if (Allscalar[i] > scalar)
					scalar = Allscalar[i];
			}
		}

		void Sent_receive_Initial_streamlines_from0(
			std::vector<Streamline>& S,
			int my_rank, boost::mpi::communicator& world) {
			world.barrier();
			int n_proc = world.size();
			if (n_proc == 1)
				return;

			std::vector<std::vector<double> > px(n_proc);
			std::vector<std::vector<double> > py(n_proc);
			std::vector<std::vector<double> > pz(n_proc);
			std::vector<std::vector<int> > E_id(n_proc);
			std::vector<std::vector<int> > S_id(n_proc);
            std::vector<std::vector<double> > pt(n_proc);
			//std::vector<std::vector<int> > proc_id(n_proc);

			// copy the data
			for (unsigned int i = 0; i < S.size(); ++i) {
				px[my_rank].push_back(S[i].getLastParticle().getP().x);
				py[my_rank].push_back(S[i].getLastParticle().getP().y);
				pz[my_rank].push_back(S[i].getLastParticle().getP().z);
				E_id[my_rank].push_back(S[i].getEid());
				S_id[my_rank].push_back(S[i].getSid());
				pt[my_rank].push_back(S[i].getLastParticle().getTime());
			}
			world.barrier();
			// Send everything to every processor
			std::vector<int> data_per_proc;
			
			Send_receive_size(static_cast<int>(px[my_rank].size()), n_proc, data_per_proc, world);
			//DEBUG::DataPerProc<int>(data_per_proc, my_rank, "After");
			Send_receive_data<double>(px, data_per_proc, my_rank, world, MPI_DOUBLE);
			Send_receive_data<double>(py, data_per_proc, my_rank, world, MPI_DOUBLE);
			Send_receive_data<double>(pz, data_per_proc, my_rank, world, MPI_DOUBLE);

			Send_receive_data<int>(E_id, data_per_proc, my_rank, world, MPI_INT);
			Send_receive_data<int>(S_id, data_per_proc, my_rank, world, MPI_INT);
            Send_receive_data<double>(pt, data_per_proc, my_rank, world, MPI_DOUBLE);
			//DEBUG::PrintVectorData<int>(S_id, my_rank, "E_id");
			S.clear();

			for (unsigned int i = 0; i < px.size(); ++i) 
				for (unsigned int j = 0; j < px[i].size(); ++j) 
					S.push_back(Streamline(E_id[i][j], S_id[i][j], Particle(vec3(px[i][j], py[i][j], pz[i][j]), pt[i][j])));
			world.barrier();
		 }

		void Send_receive_streamlines(std::vector<Streamline>& Ssend, 
			std::vector<Streamline>& Srecv, 
			boostPolygon& thisdomain,
			boost::mpi::communicator& world) {
			int my_rank = world.rank();
			int n_proc = world.size();
			Srecv.clear();

			std::vector<std::vector<double> > px(n_proc);
			std::vector<std::vector<double> > py(n_proc);
			std::vector<std::vector<double> > pz(n_proc);
			std::vector<std::vector<double> > vx(n_proc);
			std::vector<std::vector<double> > vy(n_proc);
			std::vector<std::vector<double> > vz(n_proc);
            std::vector<std::vector<double> > pt(n_proc);
			std::vector<std::vector<int> > E_id(n_proc);
			std::vector<std::vector<int> > S_id(n_proc);

			//std::vector<std::vector<int> > proc_id(n_proc);
			std::vector<std::vector<int> > p_id(n_proc);
			std::vector<std::vector<int> > Nstuck(n_proc);
			std::vector<std::vector<double> > BBlx(n_proc);
			std::vector<std::vector<double> > BBly(n_proc);
			std::vector<std::vector<double> > BBlz(n_proc);
			std::vector<std::vector<double> > BBux(n_proc);
			std::vector<std::vector<double> > BBuy(n_proc);
			std::vector<std::vector<double> > BBuz(n_proc);
			std::vector<std::vector<double> > age(n_proc);
			world.barrier();
			for (unsigned int i = 0; i < Ssend.size(); ++i) {
				px[my_rank].push_back(Ssend[i].getLastParticle().getP().x);
				py[my_rank].push_back(Ssend[i].getLastParticle().getP().y);
				pz[my_rank].push_back(Ssend[i].getLastParticle().getP().z);
				vx[my_rank].push_back(Ssend[i].getLastParticle().getV().x);
				vy[my_rank].push_back(Ssend[i].getLastParticle().getV().y);
				vz[my_rank].push_back(Ssend[i].getLastParticle().getV().z);
                pt[my_rank].push_back(Ssend[i].getLastParticle().getTime());
				E_id[my_rank].push_back(Ssend[i].getEid());
				S_id[my_rank].push_back(Ssend[i].getSid());
				//proc_id[my_rank].push_back(Ssend[i].getLastParticle().getProc());
				p_id[my_rank].push_back(Ssend[i].getLastParticle().getPid());
				Nstuck[my_rank].push_back(Ssend[i].StuckIter());
				BBlx[my_rank].push_back(Ssend[i].getBBlow().x);
				BBly[my_rank].push_back(Ssend[i].getBBlow().y);
				BBlz[my_rank].push_back(Ssend[i].getBBlow().z);
				BBux[my_rank].push_back(Ssend[i].getBBupp().x);
				BBuy[my_rank].push_back(Ssend[i].getBBupp().y);
				BBuz[my_rank].push_back(Ssend[i].getBBupp().z);
				age[my_rank].push_back(Ssend[i].getAge());
				
				
				//if (my_rank == 1) {
				//	std::cout << E_id[my_rank][i] << ", " << S_id[my_rank][i] << ", " << p_id[my_rank][i] << ", " << proc_id[my_rank][i] << ", "
				//		<< px[my_rank][i] << ", " << py[my_rank][i] << ", " << pz[my_rank][i] << ", "
				//		<< vx[my_rank][i] << ", " << vy[my_rank][i] << ", " << vz[my_rank][i] << std::endl;
				//}
				
			}
			world.barrier();
			std::vector<int> data_per_proc;
			MPI::Send_receive_size(px[my_rank].size(), n_proc, data_per_proc, world);
			//if (my_rank == 1) {
			//	for (int i = 0; i < data_per_proc.size(); ++i)
			//		std::cout << "data_per_proc " << data_per_proc[i] << std::endl;
			//}

			MPI::Send_receive_data<double>(px, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(py, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(pz, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(vx, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(vy, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(vz, data_per_proc, my_rank, world, MPI_DOUBLE);
            MPI::Send_receive_data<double>(pt, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<int>(E_id, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<int>(S_id, data_per_proc, my_rank, world, MPI_INT);
			//MPI::Send_receive_data<int>(proc_id, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<int>(p_id, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<int>(Nstuck, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<double>(BBlx, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBly, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBlz, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBux, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBuy, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBuz, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(age, data_per_proc, my_rank, world, MPI_DOUBLE);

			/*
			if (my_rank == 0) {
				for (int i = 0; i < px.size(); ++i) {
					for (int j = 0; j < px[i].size(); ++j) {
						std::cout << E_id[i][j] << ":" << S_id[i][j] << ":" << p_id[i][j] << ":" << proc_id[i][j] << std::endl;
						std::cout << "px[" << i << "][" << j << "]=" << px[i][j] <<  std::endl;
						std::cout << "py[" << i << "][" << j << "]=" << py[i][j] << std::endl;
						std::cout << "pz[" << i << "][" << j << "]=" << pz[i][j] << std::endl;

					}
				}
			}
			*/
			
			// Now loop through the data pick the ones that the other processors found 
			// that they should go on the current processor
			world.barrier();
			for (unsigned int i = 0; i < px.size(); ++i) {
				if (i == my_rank)
					continue;
				for (unsigned int j = 0; j < px[i].size(); ++j) {
					if (boost::geometry::within(boostPoint(px[i][j], py[i][j]), thisdomain )) {
						vec3 p = vec3(px[i][j], py[i][j], pz[i][j]);
						vec3 v = vec3(vx[i][j], vy[i][j], vz[i][j]);
						vec3 bl = vec3(BBlx[i][j], BBly[i][j], BBlz[i][j]);
						vec3 bu = vec3(BBux[i][j], BBuy[i][j], BBuz[i][j]);
						Streamline Stmp = Streamline(E_id[i][j], S_id[i][j],
							Particle(p, v, p_id[i][j] /*, proc_id[i][j]*/, pt[i][j]),
							bl, bu, Nstuck[i][j], age[i][j]);
						Srecv.push_back(Stmp);
					}
				}
			}
			world.barrier();
			//std::cout << my_rank << " Srecv size: " << Srecv.size() << std::endl;
		}
	 } 
}
