#pragma once
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// Includes for interpolation
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

// Includes for dD Search
//#include <CGAL/Epick_d.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// 2D search 
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_set_2.h>



#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <iostream>
#include <fstream>
#include <vector>


namespace ICHNOS {

    const double sqrt2 = std::sqrt(2.0);
	/**
	 * @brief This is a structure to hold parameters for the distribution of 
	 * particles around the wells
	 * 
	 */
	struct WellOptions {
		//! The total number of particles around the wells 
		int Nparticles;
		
		//! The number of layers that the particles are orginized.
		int Nlayer;
		
		//! The particles are distributed around some distance Radius
		double Radius;
	};


	/**
	 * @brief A structure that holds a vector and some operators for the vector
	 * 
	 */
	struct vec3
	{
		double x = 0; //! The x coordinate
		double y = 0; //! The y coordinate
		double z = 0; //! The z coordinate

		//! The default constructor sets everything to zero
		vec3(double x = 0, double y = 0, double z = 0)
			:
			x(x), y(y), z(z)
		{}

		void operator=(const vec3& a) {
		    x = a.x;
		    y = a.y;
		    z = a.z;
		}

        void operator=(const double& a) {
		    x = a;
		    y = a;
		    z = a;
        }

		//! Element wise vector multiplication operator
		vec3 operator*(const vec3& a) const {
			return vec3(a.x * x, a.y * y, a.z * z);
		}

		//! Multiplication with a scalar. The scale should go to the right side 
		vec3 operator*(const double& a) const {
			return vec3(a * x, a * y, a * z);
		}

		//! add operator
		vec3 operator+(const vec3& a) const {
			return vec3(a.x + x, a.y + y, a.z + z);
		}

		//! Subtraction operator. Subtracts the other vector from this
		vec3 operator-(const vec3& a) const {
			//return vec3(a.x - x, a.y - y, a.z - z);
			return vec3(x - a.x, y - a.y, z - a.z);
		}

		//! Calculates the length of the vector
		double len() {
			return std::sqrt(x * x + y * y + z * z);
		}

		//! Checks whether the vector is almost zero
		bool isZero() {
			bool tf = len() < 0.0000000001;
			return tf;
		}

		//! Some methods sets the vector to -99999 when the calculations failed. This methods checks is the vector has this value 
		bool isInvalid() {
			bool tf = std::sqrt((x + 99999) * (x + 99999) + (y + 99999) * (y + 99999) + (z + 99999) * (z + 99999)) < 0.0000000001;
			return tf;
		}

		//! make the length of the vector 1.
		vec3 normalize() {
			double length = len();
			return (vec3(x / length, y / length, z / length));
		}

		//! reset the values of the vector back to 0
		void zero() {
			x = 0.0;
			y = 0.0;
			z = 0.0;
		}

		double distance(double px, double py, double pz) {
			return std::sqrt((x - px) * (x - px) + (y - py) * (y - py) + (z - pz) * (z - pz));
		}

		double angle(vec3& b){
			double ab = x*b.x + y*b.y + z*b.z;
			double abs_a = this->len();
			double abs_b = b.len();
			return std::acos(ab/(abs_a*abs_b));
		}
	};

    struct TimeData{
        double tm;
        double tm_tmp;
        int idx1;
        int idx2;
        double t;
    };

    struct helpVars{
        vec3 pp;
        vec3 vv;
        vec3 ll;
        vec3 uu;
        double adaptStepSize;
        double actualStep;
        double actualStepTime;
        double diameter;
        double ratio;
        bool bIsInitialized;
        TimeData td;
    };


	struct NPSAT_data {
		int proc = -9;
		int id = -9;
		double diameter = 0;
		double ratio = 0;
		vec3 v;
	};

    struct StepOptions {
        /// This is the step size with units of length. For RK45 this is the initial step size
        double StepSize;
        /// This is the step size with units of Time. For RK45 this is the initial step size
        double StepSizeTime;
        /// This is the number of steps to divide the size of an element
        double nSteps;
        /// This is the number of steps to take within a time step
        double nStepsTime;
        /// This is the minimum step size at the exit side of the particles
        double minExitStepSize;
        /// The direction of particle tracking
        double dir;
    };

	struct elev_data{
		double top = 0;
		double bot = 0;
		int id = -9;
	};

	// 3D Triangulation
	typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
	typedef CGAL::Triangulation_vertex_base_with_info_3<NPSAT_data, K>  Vb;
	typedef CGAL::Delaunay_triangulation_cell_base_3<K>	Cb;
	typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
	typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>	cgal_Delaunay;
	typedef cgal_Delaunay::Point cgal_point;
	typedef cgal_Delaunay::Vertex_handle vertex_handle;
	// typedef cgal_Delaunay::Cell_handle cell_handle;

	// Interpolation
	// typedef CGAL::Delaunay_triangulation_2<K> cgal_Delaunay_2;
	// typedef CGAL::Interpolation_traits_2<K>	interp_traits;
	// typedef K::FT coord_type;
	//

	// typedef std::map<cgal_point_2, coord_type, K::Less_xy_2> coord_map;
	// typedef CGAL::Data_access<coord_map>	value_access;

	// dD Spatial Searching
	typedef K::Point_3	cgal_point_3;
	typedef boost::tuple<cgal_point_3, int> pnt_int;
	typedef CGAL::Search_traits_3<K> Traits_base;
	typedef CGAL::Search_traits_adapter<pnt_int, CGAL::Nth_of_tuple_property_map<0,pnt_int>, Traits_base> search_traits;
	typedef CGAL::Fuzzy_iso_box<search_traits> Fuzzy_iso_box;
	typedef CGAL::Kd_tree<search_traits> search_tree_int;
    typedef CGAL::Orthogonal_k_neighbor_search<search_traits> K_neighbor_search_int;

	

	// 2D searching
	typedef CGAL::Triangulation_vertex_base_with_info_2<elev_data, K>   Vb2D;
	typedef CGAL::Triangulation_data_structure_2<Vb2D>  Tds2D;
	typedef CGAL::Point_set_2<K,Tds2D>::Vertex_handle Vertex_handle2D;
	typedef CGAL::Point_set_2<K,Tds2D>  PointSet2;
	typedef K::Point_2 cgal_point_2;


	// Nstates Nmonths Nvelocities
	typedef std::vector<std::vector<std::vector<vec3>>> stoch_vel;

	class Stochastic_Velocity {
	public:
		Stochastic_Velocity();
		Stochastic_Velocity(int nstates, int nperiods);
		void initStatesPeriods(int nstates, int nperiods);
		void initStatesPeriods();
		void addValue(int istate, int iperiod, vec3 value);
		int nValues(int istate, int iperiod);
		void getValue(int istate, int iperiod, int r, vec3& value);

	private:
		int Nstates;
		int Nperiods;
		stoch_vel VelocityField;
	};

	Stochastic_Velocity::Stochastic_Velocity()
		:
		Nstates(-9),
		Nperiods(-9)
	{}

	Stochastic_Velocity::Stochastic_Velocity(int nstates, int nperiods)
		:
		Nstates(nstates),
		Nperiods(nperiods)
	{}

	void Stochastic_Velocity::initStatesPeriods() {
		VelocityField.clear();
		std::vector< std::vector<vec3> > mon_vec;
		std::vector<vec3> vel_vec;
		for (int i = 0; i < Nperiods; i++) {
			mon_vec.push_back(vel_vec);
		}
		for (int i = 0; i < Nstates; i++) {
			VelocityField.push_back(mon_vec);
		}
	}

	void Stochastic_Velocity::initStatesPeriods(int nstates, int nperiods) {
		Nstates = nstates;
		Nperiods = nperiods;
		initStatesPeriods();
	}

	void Stochastic_Velocity::addValue(int istate, int iperiod, vec3 value) {
		VelocityField[istate][iperiod].push_back(value);
	}
	int Stochastic_Velocity::nValues(int istate, int iperiod) {
		return static_cast<int>(VelocityField[istate][iperiod].size());
	}

	void Stochastic_Velocity::getValue(int istate, int iperiod, int r, vec3& value) {
		value = VelocityField[istate][iperiod][r];
	}



	// Data structure for stochastic simulation
	struct STOCH_data {
		int proc = -9;
		int id = -9;
		double diameter = 0;
		double ratio = 0;
		Stochastic_Velocity v;
	};

	// Stochastic tree
	typedef boost::tuple<cgal_point_3, STOCH_data> pnt_stoch;
	typedef CGAL::Search_traits_adapter<pnt_stoch, CGAL::Nth_of_tuple_property_map<0, pnt_stoch>, Traits_base> search_traits_stoch;
	typedef CGAL::Fuzzy_iso_box<search_traits_stoch> Fuzzy_iso_box_stoch;
	typedef CGAL::Kd_tree<search_traits_stoch> search_tree_stoch;

    enum class infoType {Nelem, Nface};
    enum class coordDim {vx, vy, vz};
    enum class TimeInterpType {NEAREST, LINEAR};
    enum class MeshVelInterpType{ELEMENT, NODE, FACE, UNKNOWN};



    struct Pnt_info{
        int proc = -9;
        int id = -9;
        double diameter = 0;
        double ratio = 0;
    };

    // Transient tree
    typedef boost::tuple<cgal_point_3, Pnt_info> pnt_with_info;
    typedef CGAL::Search_traits_adapter<pnt_with_info, CGAL::Nth_of_tuple_property_map<0, pnt_with_info>, Traits_base> search_traits_pnt_info;
    typedef CGAL::Fuzzy_iso_box<search_traits_pnt_info> Fuzzy_iso_box_info;
    typedef CGAL::Kd_tree<search_traits_pnt_info> search_tree_info;
    typedef CGAL::Orthogonal_k_neighbor_search<search_traits_pnt_info> K_neighbor_search;


    struct Pnt_IWFM_info{
        int proc = -9;
        int tri_id = -9;
        int elem_id = -9;
        double diameter = 0;
        std::vector<double> height;
    };
    // IWFM tree
    typedef boost::tuple<cgal_point_3, Pnt_IWFM_info> pnt_iwfm;
    typedef CGAL::Search_traits_adapter<pnt_iwfm, CGAL::Nth_of_tuple_property_map<0, pnt_iwfm>, Traits_base> search_traits_pnt_iwfm;
    typedef CGAL::Fuzzy_sphere<search_traits_pnt_iwfm> Fuzzy_sphere_iwfm;
    typedef CGAL::Kd_tree<search_traits_pnt_iwfm> search_tree_iwfm;

    class VelTR{
    public:
        VelTR();
        void init(int np, int nt, int dim_in = 3);
        void init(int np, int dim_in = 3);
        void setVELvalue(double v, int pnt, int step, coordDim dim);
        void setTSvalue(double v, int step);
        void setTSvalue(std::vector<double>& TS_in);
        void findIIT(double x, int &i1, int &i2, double &t, double &x_tmp);
        vec3 getVelocity(int pnt, int i1, int i2, double t);
        void getVelocity(std::vector<int> &pnts, int i1, int i2, double t, std::vector<vec3> &vel_out);
        double getTSvalue(int idx);
        void setNrepeatDays(double n){nDaysRepeat = n;}
        void setTimeInterpolationType(TimeInterpType tip_in){tip = tip_in;}
        void setnSteps(int n){nSteps = n;}
        int getNpoints(){return nPoints;}
        double stepTimeUpdate(helpVars &pvlu, StepOptions &stepOpt);
    private:
        std::vector<std::vector<double>> VX;
        std::vector<std::vector<double>> VY;
        std::vector<std::vector<double>> VZ;

        std::vector<double> TS;

        int nPoints;
        int nSteps;
        int nDims;
        void findTimeStepIndex(int &i, int &ii, double x);
        void checkIndices(int &i1, int &i2);

        double nDaysRepeat = 0;
        TimeInterpType tip;
        bool bIsInitialzed = false;

    };

    VelTR::VelTR()
        :
        nPoints(-9),
        nSteps(-9)
    {}

    void VelTR::init(int np, int nt, int dim_in) {
        if (bIsInitialzed)
            return;
        nPoints = np;
        nSteps = nt;
        nDims = dim_in;
        VX.clear();
        VY.clear();
        VZ.clear();
        VX.resize(nPoints, std::vector<double>(nSteps, 0));
        if (nDims >= 2)
            VY.resize(nPoints, std::vector<double>(nSteps, 0));
        if (nDims >= 3)
            VZ.resize(nPoints, std::vector<double>(nSteps, 0));
        //TS.clear();
        //TS.resize(nSteps, 0);
        bIsInitialzed = true;
    }

    void VelTR::init(int np, int dim_in){
        init(np, nSteps, dim_in);
    }

    void VelTR::setVELvalue(double v, int pnt, int step, coordDim dim){
        if (pnt >= nPoints || step >= nSteps || pnt < 0  || step < 0)
            return;

        if (dim == coordDim::vx && nDims >= 1){
            VX[pnt][step] = v;
        }
        else if (dim == coordDim::vy && nDims >= 2){
            VY[pnt][step] = v;
        }
        else if (dim == coordDim::vz && nDims >= 1){
            VZ[pnt][step] = v;
        }
    }

    void VelTR::setTSvalue(double v, int step){
        if (step >= nSteps || step < 0)
            return;
        TS[step] = v;
    }
    void VelTR::setTSvalue(std::vector<double> &TS_in) {
        TS = TS_in;
    }

    void VelTR::findTimeStepIndex(int &i, int &ii, double x){
        int mid = (i + ii)/2;
        if (x <= TS[mid])
            ii = mid;
        else
            i = mid;

        if (ii - i <= 1)
            return;
        else
            findTimeStepIndex(i, ii, x);
    }

    void VelTR::checkIndices(int &i1, int &i2) {
        if (i1 >= TS.size()-1){
            i1 = TS.size()-2;
        }
        if (i2 >= TS.size()-1){
            i2 = TS.size()-2;
        }
    }

    void VelTR::findIIT(double x, int &i1, int &i2, double &t, double &x_tmp){
        x_tmp = x;
        if (x <= TS[0]){
            if (nDaysRepeat == 0){
                i1 = 0;
                i2 = 0;
                t = 0.0;
                return;
            }
            else{
                while (x_tmp <= TS[0]){
                    x_tmp += nDaysRepeat;
                }
            }
        }
        if (x >= TS[TS.size()-1]){
            if (nDaysRepeat == 0){
                i1 = TS.size()-2;
                i2 = TS.size()-2;
                t = 1.0;
                return;
            }
            else{
                while (x_tmp >= TS[TS.size()-1]){
                    x_tmp -= nDaysRepeat;
                }
            }

        }
        i1 = 0;
        i2 = TS.size()-1;
        findTimeStepIndex(i1, i2, x_tmp);
        t = (x_tmp - TS[i1]) / (TS[i2] - TS[i1]);
        checkIndices(i1, i2);
    }

    vec3 VelTR::getVelocity(int pnt, int i1, int i2, double t){
        if (pnt >= nPoints){
            std::cout << "The requested id:" << pnt << " exceeds the number of velocity points: " << nPoints << std::endl;
            std::cout << "the returned velocity will be invalid" << std::endl;
            return vec3(-99999,-99999,-99999);
        }
        if (tip == TimeInterpType::LINEAR){
            if (nDims == 3){
                vec3 v1(VX[pnt][i1], VY[pnt][i1], VZ[pnt][i1]);
                vec3 v2(VX[pnt][i2], VY[pnt][i2], VZ[pnt][i2]);
                return v1*(1-t) + v2*t;
            }
            else if (nDims == 2){
                vec3 v1(VX[pnt][i1], VY[pnt][i1], 0.0);
                vec3 v2(VX[pnt][i2], VY[pnt][i2], 0.0);
                return v1*(1-t) + v2*t;
            }
            else if (nDims == 1){
                vec3 v1(VX[pnt][i1], 0.0, 0.0);
                vec3 v2(VX[pnt][i2], 0.0, 0.0);
                return v1*(1-t) + v2*t;
            }
        }
        else{
            if (nDims == 3)
                return vec3(VX[pnt][i1], VY[pnt][i1], VZ[pnt][i1]);
            else if (nDims == 2)
                return vec3(VX[pnt][i1], VY[pnt][i1], 0.0);
            else if (nDims == 1)
                return vec3(VX[pnt][i1], 0.0, 0.0);
        }
        return vec3(-99999,-99999,-99999);
    }

    void VelTR::getVelocity(std::vector<int> &pnts, int i1, int i2, double t, std::vector<vec3> &vel_out) {
        vel_out.clear();
        for (unsigned int i = 0; i <pnts.size(); ++i){
            vel_out.push_back(getVelocity(pnts[i],i1,i2,t));
        }
    }

    double VelTR::getTSvalue(int idx) {
        if (idx >= 0 && idx < static_cast<int>(TS.size()))
            return TS[idx];
        else
            return -9.0;
    }

    double VelTR::stepTimeUpdate(helpVars &pvlu, StepOptions &stepOpt){
        if (nSteps == 0){
            return 999999999;
        }
        double stepTime, tm_1, tm_2;
        double dt = 1/stepOpt.nStepsTime;
        if (pvlu.td.idx1 == pvlu.td.idx2){
            if (stepOpt.dir > 0){
                tm_1 = getTSvalue(pvlu.td.idx1);
                tm_2 = getTSvalue(pvlu.td.idx1+1);
                if (tm_2 < 0){
                    tm_2 = tm_1;
                    tm_1 = getTSvalue(pvlu.td.idx1-1);
                }
            }
            else{
                tm_1 = getTSvalue(pvlu.td.idx1-1);
                tm_2 = getTSvalue(pvlu.td.idx1);
                if (tm_1 < 0){
                    tm_1 = tm_2;
                    tm_2 = getTSvalue(pvlu.td.idx1+1);
                }
            }
        }
        else{
            tm_1 = getTSvalue(pvlu.td.idx1);
            tm_2 = getTSvalue(pvlu.td.idx2);
        }

        double tmp_step = dt*(tm_2 - tm_1);
        stepTime = pvlu.vv.len() * tmp_step;
        /*double end_time = pvlu.td.tm + stepOpt.dir*tmp_step;
        if (stepOpt.dir > 0){
            if (std::abs(pvlu.td.tm - tm_2) < 0.25*tmp_step){
                if (pvlu.td.idx2 + 1 < nSteps){
                    tm_2 = getTSvalue(pvlu.td.idx2+1);
                }
            }
            if (end_time > tm_2 && pvlu.td.idx2 < nSteps - 1){
                stepTime = pvlu.vv.len() * (tm_2 - pvlu.td.tm);
            }
            else{
                stepTime = pvlu.vv.len() * tmp_step;
            }
        }
        else{
            if (end_time < tm_1 && pvlu.td.idx1 > 0){
                stepTime = pvlu.vv.len() * (pvlu.td.tm - tm_1);
            }
            else{
                stepTime = pvlu.vv.len() * tmp_step;
            }
        }*/

        return stepTime;
    }


	enum class interpType { INGORE, SCALAR, CLOUD};

	enum class SolutionMethods
	{
		Euler,
		RK2,
		RK4,
		RK45,
		INVALID
	};

	

	SolutionMethods castMethod2Enum(std::string method) {
		std::map < std::string, SolutionMethods> SolutionMethodMap;
		std::map < std::string, SolutionMethods>::iterator it;
		SolutionMethodMap.insert(std::pair<std::string, SolutionMethods>("Euler", SolutionMethods::Euler));
		SolutionMethodMap.insert(std::pair<std::string, SolutionMethods>("RK2", SolutionMethods::RK2));
		SolutionMethodMap.insert(std::pair<std::string, SolutionMethods>("RK4", SolutionMethods::RK4));
		SolutionMethodMap.insert(std::pair<std::string, SolutionMethods>("RK45", SolutionMethods::RK45));
		it = SolutionMethodMap.find(method);
		if (it != SolutionMethodMap.end())
			return it->second;
		else {
			return SolutionMethods::INVALID;
		}
	}

	std::string castMethod2String(SolutionMethods m) {
		std::map <SolutionMethods, std::string> SolutionMethodMap;
		std::map <SolutionMethods, std::string>::iterator it;
		SolutionMethodMap.insert(std::pair<SolutionMethods, std::string>(SolutionMethods::Euler, "Euler"));
		SolutionMethodMap.insert(std::pair<SolutionMethods, std::string>(SolutionMethods::RK2, "RK2"));
		SolutionMethodMap.insert(std::pair<SolutionMethods, std::string>(SolutionMethods::RK4, "RK4"));
		SolutionMethodMap.insert(std::pair<SolutionMethods, std::string>(SolutionMethods::RK45, "RK45"));
		it = SolutionMethodMap.find(m);
		if (it != SolutionMethodMap.end())
			return it->second;
		else {
			return "INVALID";
		}
	}

	enum class ExitReason {
		CHANGE_PROCESSOR,
		EXIT_TOP,
		EXIT_BOTTOM,
		EXIT_SIDE,
		INIT_OUT,
		STUCK,
		NO_EXIT,
		MAX_INNER_ITER,
		FIRST_POINT_GHOST,
		FAR_AWAY,
		MAX_AGE,
		ATTRACT,
        NOT_IN_SUBDOMAIN,
		NO_REASON
	};

	std::string castExitReasons2String(ExitReason er) {
		std::map <ExitReason, std::string> ExitReasonsMap;
		std::map <ExitReason, std::string>::iterator it;
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::CHANGE_PROCESSOR, "CHANGE_PROCESSOR"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::EXIT_TOP, "EXIT_TOP"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::EXIT_BOTTOM, "EXIT_BOTTOM"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::EXIT_SIDE, "EXIT_SIDE"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::INIT_OUT, "INIT_OUT"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::STUCK, "STUCK"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::NO_EXIT, "NO_EXIT"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::MAX_INNER_ITER, "MAX_INNER_ITER"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::FIRST_POINT_GHOST, "FIRST_POINT_GHOST"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::FAR_AWAY, "FAR_AWAY"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::MAX_AGE, "MAX_AGE"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::ATTRACT, "ATTRACT"));
        ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::NOT_IN_SUBDOMAIN, "NOT_IN_SUBDOMAIN"));
		ExitReasonsMap.insert(std::pair<ExitReason, std::string>(ExitReason::NO_REASON, "NO_REASON"));
		it = ExitReasonsMap.find(er);
		if (it != ExitReasonsMap.end())
			return it->second;
		else {
			return "INVALID";
		}
	}

	ExitReason castExitReasons2Enum(std::string exitreason) {
		std::map < std::string, ExitReason> ExitReasonMap;
		std::map < std::string, ExitReason>::iterator it;
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("CHANGE_PROCESSOR", ExitReason::CHANGE_PROCESSOR));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("EXIT_BOTTOM", ExitReason::EXIT_BOTTOM));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("EXIT_SIDE", ExitReason::EXIT_SIDE));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("EXIT_TOP", ExitReason::EXIT_TOP));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("FAR_AWAY", ExitReason::FAR_AWAY));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("MAX_AGE", ExitReason::MAX_AGE));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("ATTRACT", ExitReason::ATTRACT));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("FIRST_POINT_GHOST", ExitReason::FIRST_POINT_GHOST));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("INIT_OUT", ExitReason::INIT_OUT));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("MAX_INNER_ITER", ExitReason::MAX_INNER_ITER));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("NO_EXIT", ExitReason::NO_EXIT));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("NO_REASON", ExitReason::NO_REASON));
		ExitReasonMap.insert(std::pair<std::string, ExitReason>("STUCK", ExitReason::STUCK));
        ExitReasonMap.insert(std::pair<std::string, ExitReason>("NOT_IN_SUBDOMAIN", ExitReason::NOT_IN_SUBDOMAIN));
		it = ExitReasonMap.find(exitreason);
		if (it != ExitReasonMap.end())
			return it->second;
		else {
			return ExitReason::NO_REASON;
		}
	}



	struct AdaptStepOptions{
        /// When the method is adaptive this is the maximum step size
        double MaxStepSize;
        /// When the method is adaptive this is the minimum step size
        double MinStepSize;
        /// When the method is adaptive this is the tolerance
        double ToleranceStepSize;
        /// When the method is adaptive this modifies the rate of change of the step size.
        double increaseRateChange;
        /// When the method is adaptive this limits the upper threshold of step decrease.
        /// For example if the decrease is lower that 1 but very close to one then it may be possible that
        /// the algorithm will repeat the decreasing of the step size many times by a very small amount.
        double limitUpperDecreaseStep;
	};


	struct ParticleOptions {
		
		SolutionMethods method;
        StepOptions StepOpt;
        AdaptStepOptions AdaptOpt;
		// The maximum number of iteratiosn that each streamline is alloweded not to expand
		int StuckIterations;

		// If update step size is set to 1 then the step size is determined 
		// by the extend of the bounding box that the particle currently is.
		// If its zero the the StepSize is used unless the method is RK45
		int UpdateStepSize;

		// The direction of the particle tracing
		double Direction;
		// The maximum steps for each streamline
		int MaxIterationsPerStreamline;
		// The maximum number of exchanges between processors
		int MaxProcessorExchanges;

		double AgeLimit;

		std::string ParticleFile;
		std::string WellFile;
		int ParticlesInParallel;
		std::string OutputFile;
		bool printH5;
		bool printASCII;
		// During gathering puut all streamlines into one file.
		// If false each file will contain ParticlesInParallel streamlines per file
		bool gatherOneFile;
		std::string configfile;
		// Number of realizations
		int Nrealizations;
		int Nthreads;
		bool RunAsThread;
        double OutputFrequency;
        int nLeadingZeros;

		bool bIsTransient = false;
	};

	struct DomainOptions {
		std::string polygonFile;
		std::string TopElevationFile;
		std::string BottomeElevationFile;
		std::string AttractorsFile;
		double TopRadius;
		double BotRadius;
		double TopPower;
		double BotPower;
		double AttractRadius;
		std::string processorDomainFile;
		//std::string expandedDomainFile;
		int myRank;
		int nProc;
		bool RunAsThread;
	};


	class Particle {
	public:
		Particle(vec3 p, double rt = 0);
		Particle(vec3 p, vec3 v, int pid/*, int proc*/, double rt = 0);
		Particle(vec3 p, vec3 v, /*int proc,*/  Particle prev, double rt = 0);
		vec3 getP() const { return p; }
		vec3 getV() const { return v; }
		double getTime() const {return ptime; }
		//double getAge() const { return age; }
		//void SetAge(double a) { age = a; }
		void SetV(vec3 vel) { v = vel; }
		int getPid() const { return pid; }
		void displayAsVEX(bool printAttrib);
		//void SetProc(int procId) { proc = procId; }
		//int getProc() { return proc; }
	private:
		int pid;
		vec3 p;
		vec3 v;
		double ptime;
		//double age;
		//int proc;
	};

	Particle::Particle(vec3 p, double rt)
		:
            p(p),
            ptime(rt)
	{
		pid = 0;
		//age = 0;
		//proc = -9;
	}

	Particle::Particle(vec3 p, vec3 v, int pid/*, int proc*/, double rt)
		:
		pid(pid),
		p(p),
		v(v),
		ptime(rt)
		//proc(proc)
	{
		//age = 0;
	}

	Particle::Particle(vec3 p, vec3 v, /*int proc,*/ Particle prev, double rt)
		:
		p(p),
		v(v),
        ptime(rt)
		//proc(proc)
	{
		pid = prev.getPid() + 1;
		//double d = (p - prev.getP()).len();
		//age = prev.getAge() + d / prev.getV().len();
	}

	void Particle::displayAsVEX(bool printAttrib) {
		vec3 vn = v * (1 / v.len());
		std::cout << std::setprecision(3) << std::fixed
		          << "p = addpoint(0,{" << p.x << "," << p.z << "," << p.y << "});";
		if (printAttrib)
			std::cout
                    << std::setprecision(6) << std::fixed << " setpointattrib(0,'N',p,{" << vn.x << "," << vn.z << "," << vn.y << "},'set');";
		std::cout << std::endl;
	}


	class Streamline {
	public:
		Streamline(int Eid_in, int Sid_in, Particle particle);
		Streamline(int Eid_in, int Sid_in, Particle particle, vec3 BL, vec3 BU, int N, double age);
		void InitialVelocity(vec3 v) {
			SL.back().SetV(v);
		}
		void AddParticle(Particle p);
		Particle getLastParticle() { return SL.back(); }
		Particle getParticleBeforeLast(){
			if (SL.size() > 1)
				return SL[SL.size()-2];
			else
			{
				return Particle(vec3());
			}
		}
		void Close(ExitReason ER);
		int StuckIter() { return IterationsWithoutExpansion; }
		int getEid() { return Eid; }
		int getSid() { return Sid; }
		vec3 getBBlow() { return BL; }
		vec3 getBBupp() { return BU; }
		unsigned int size() { return SL.size(); }
		Particle getParticle(int i) { return SL.at(i); }
		double getAge() { return age; }
        helpVars PVLU;
        ExitReason getExitReason(){return exitreason;}
        bool printIt(){return PrintThis;}

	private:
		int Eid;
		int Sid;
		std::vector<Particle > SL;
		ExitReason exitreason = ExitReason::NO_REASON;
        bool PrintThis = false;
		vec3 BL;
		vec3 BU;
		int IterationsWithoutExpansion = 0;
		// This is the time the particle has spend in the velocity field
		double age = 0;
	};

	Streamline::Streamline(int Eid_in, int Sid_in, Particle particle)
		:
		Eid(Eid_in),
		Sid(Sid_in)
	{
		SL.clear();
		SL.push_back(particle);
		BL = particle.getP();
		BU = particle.getP();
		IterationsWithoutExpansion = 0;
		age = 0;
	}

	Streamline::Streamline(int Eid_in, int Sid_in, Particle particle, vec3 BL, vec3 BU, int N, double age_in)
		:
		Eid(Eid_in),
		Sid(Sid_in),
		BL(BL),
		BU(BU),
		IterationsWithoutExpansion(N),
		age(age_in)
	{
		SL.clear();
		SL.push_back(particle);
	}

	void Streamline::AddParticle(Particle p) {
		vec3 p_prev = getLastParticle().getP();
		double vv = getLastParticle().getV().len();
		SL.push_back(p);
		vec3 pp = p.getP();
		double dst = (pp - p_prev).len();
		age += dst / vv;

		bool expand = false;
		if (pp.x < BL.x) { BL.x = pp.x; expand = true; }
		if (pp.x > BU.x) { BU.x = pp.x; expand = true; }
		if (pp.y < BL.y) { BL.y = pp.y; expand = true; }
		if (pp.y > BU.y) { BU.y = pp.y; expand = true; }
		if (pp.z < BL.z) { BL.z = pp.z; expand = true; }
		if (pp.z > BU.z) { BU.z = pp.z; expand = true; }
		if (expand)
			IterationsWithoutExpansion = 0;
		else
			IterationsWithoutExpansion++;
	}

	void Streamline::Close(ExitReason ER) {
		exitreason = ER;
        if (ER == ExitReason::NOT_IN_SUBDOMAIN)
            PrintThis = false;
        else
            PrintThis = true;
    }

	std::vector<std::vector<double> > getRK45coef() {
		std::vector<std::vector<double> > RKcoef;
		RKcoef.push_back(std::vector<double>{1.0 / 4.0, 1.0 / 4.0});
		RKcoef.push_back(std::vector<double>{3.0 / 8.0, 3.0 / 32.0, 9.0 / 32.0});
		RKcoef.push_back(std::vector<double>{12.0 / 13.0, 1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0});
		RKcoef.push_back(std::vector<double>{1.0, 439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0});
		RKcoef.push_back(std::vector<double>{1.0 / 2.0, -8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0});
		RKcoef.push_back(std::vector<double>{1.0, 25.0 / 216.0, 1408.0 /2565.0, 2197.0 / 4101.0, -1.0 / 5.0});
		RKcoef.push_back(std::vector<double>{1.0, 16.0 / 135.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 /55.0});

		return RKcoef;
	}

	// A temporary structure that is used durign the gathering phase
	struct gPart {
		vec3 p;
		vec3 v;
		double age;
	};

	struct gStream {
		std::map<int, gPart> particles;
		ExitReason ex;
	};

	//				 Eid		  Sid	
	typedef std::map<int, std::map<int, gStream > > streamlineMap;
}

//typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Adaptor<double, ICHNOS::pointCloud< ICHNOS::vec3 > >, ICHNOS::pointCloud<ICHNOS::vec3 >, 3 > nano_kd_tree_vector;
//typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Adaptor<double, ICHNOS::pointCloud< double > >, ICHNOS::pointCloud<double >, 3 > nano_kd_tree_scalar;
//typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Adaptor<double, ICHNOS::pointCloud< int > >, ICHNOS::pointCloud<int >, 3 > nano_kd_tree_int;

typedef boost::geometry::model::d2::point_xy<double> boostPoint;
typedef boost::geometry::model::polygon<boostPoint> boostPolygon;

namespace ICHNOS {

	class multiPoly {
	public:
		multiPoly() {};
		bool readfromFile(std::string filename);
		bool is_point_in(double x, double y);
	private:
		std::vector<boostPolygon> outside;
		std::vector<boostPolygon> inside;
	};


	// The format is
	// N polygons
	// Repeat the following N times
	// Np points of the polygons 1 for polygon or 0 if its a hole
	// Read Np coordinates x,y. Orientation doesn't matter and will be corrected inside the code
	bool multiPoly::readfromFile(std::string filename) {
		std::ifstream datafile(filename.c_str());
		if (!datafile.good()) {
			std::cout << "Can't open the file" << filename << std::endl;
			return false;
		}
		else {
			int Np, orien;
			double x, y;
			std::string line;
			while (getline(datafile, line))
			{
				std::istringstream inp(line.c_str());
				boostPolygon polygon;
				std::vector<boostPoint> polygonPoints;
				inp >> Np;
				inp >> orien;
				for (int j = 0; j < Np; ++j) {
					getline(datafile, line);
					std::istringstream inp1(line.c_str());
					inp1 >> x;
					inp1 >> y;
					polygonPoints.push_back(boostPoint(x, y));
				}
				boost::geometry::assign_points(polygon, polygonPoints);
				boost::geometry::correct(polygon);
				if (orien == 0) {
					outside.push_back(polygon);
				}
				else {
					inside.push_back(polygon);
				}

			}
			datafile.close();
			return true;
		}
	}

	bool multiPoly::is_point_in(double x, double y) {
		boostPoint pnt(x, y);
		for (unsigned int i = 0; i < outside.size(); ++i) {
			if (boost::geometry::within(pnt, outside[i])) {
				return false;
			}
		}
		for (unsigned int i = 0; i < inside.size(); ++i) {
			if (boost::geometry::within(pnt, inside[i])) {
				return true;
			}
		}
		return false;
	}

	class SingletonGenerator {
		static SingletonGenerator* _instance;
		//unsigned int seed;

		SingletonGenerator() {
			generetor.seed(std::time(0));
			//std::cout << "constructor" << std::endl;
            //std::cout << normalDistribution.mean() << std::endl;
            //std::cout << normalDistribution.sigma() << std::endl;
		}

		boost::random::uniform_real_distribution<float> uniformDistribution;
		boost::random::normal_distribution<float> normalDistribution;
		boost::mt19937 generetor;

	public:
		//void printSeed() {
		//	std::cout << "Seed: " << seed << std::endl;
		//}
		static SingletonGenerator* getInstance() {
			if (!_instance)
				_instance = new SingletonGenerator;
			return _instance;
		}

		float randomNumber() {
			return uniformDistribution(generetor);
		}

		float randomNumber(float min, float max) {
			return min + (max - min) * randomNumber();
		}

		int randomNumber(int min, int max) {
			return static_cast<int>(randomNumber(static_cast<float>(min), static_cast<float>(max)));
			//float r;
			//int rint;
			//while (true) {
			//	r = randomNumber(static_cast<float>(min), static_cast<float>(max));
			//	rint = static_cast<int>(r);
			//	if (rint >= min && rint <= max)
			//		break;
			//}

			//return rint;
		}
		float randomNormal(){
            //std::cout << normalDistribution.mean() << std::endl;
            //std::cout << normalDistribution.sigma() << std::endl;
            return  normalDistribution(generetor);
		}
	};

	class TransitionProbabilityMatrix {
	public:
		TransitionProbabilityMatrix() {};
		bool readData(std::string filename);
		void setNstates(int n);
		int nextState(int previousState, vec3 p, double r);
	private:
		int Nstates = 0;
		int Nmatrices = 0;
		std::vector<boostPolygon> TPMpolygons;
		std::vector<std::vector<std::vector<double>>> TPMatrices;
	};

	void TransitionProbabilityMatrix::setNstates(int n) {
		Nstates = n;
	}

	int TransitionProbabilityMatrix::nextState(int previousState, vec3 p, double r) {
		int ipoly = 0;
		int newstate = 0;
		for (int i = 0; i < Nmatrices; i++) {
			bool tf = boost::geometry::within(boostPoint(p.x, p.y), TPMpolygons[i]);
			if (tf) {
				ipoly = i;
				break;
			}
		}

		for (int i = 0; i < Nstates; i++) {
			if (r < TPMatrices[ipoly][previousState][i]){
				newstate = i;
				break;
			}
		}
		return newstate;
	}

	bool TransitionProbabilityMatrix::readData(std::string filename) {
		if (Nstates == 0) {
			std::cout << "The number ofstates is zero" << std::endl;
			return false;
		}

		std::ifstream datafile(filename.c_str());
		if (!datafile.good()) {
			std::cout << "Can't open the file " << filename << std::endl;
			return false;
		}
		else {
			std::string line;
			double x, y;
			int nverts;
			{
				getline(datafile, line);
				std::istringstream inp(line.c_str());

				inp >> Nmatrices;
			}
			for (int imat = 0; imat < Nmatrices; imat++) {
				TPMatrices.push_back(std::vector<std::vector<double>>(Nstates, std::vector<double>(Nstates)));
				{
					{
						getline(datafile, line);
						std::istringstream inp(line.c_str());
						inp >> nverts;
					}
					boostPolygon polygon;
					std::vector<boostPoint> polygonPoints;
					for (int j = 0; j < nverts; j++) {
						getline(datafile, line);
						std::istringstream inp(line.c_str());
						inp >> x;
						inp >> y;
						polygonPoints.push_back(boostPoint(x, y));
					}
					boost::geometry::assign_points(polygon, polygonPoints);
					boost::geometry::correct(polygon);
					TPMpolygons.push_back(polygon);

					for (int i = 0; i < Nstates; i++) {
						getline(datafile, line);
						std::istringstream inp(line.c_str());
						for (int j = 0; j < Nstates; j++) {
							inp >> x;
							TPMatrices[imat][i][j] = x;
						}
					}
				}
			}
			datafile.close();
			return true;
		}
	}


    class CellGraph{
    public:
        CellGraph(){};
        bool readGraphFile(std::string filename);
        bool getNearVelocities(vec3& p, std::vector<int>& velIds);

    private:
        search_tree_int Tree;
        //std::vector<vec3> XYZ;
        std::vector<std::vector<int>> neighElem;
        std::vector<std::vector<int>> velids;
    };

    bool CellGraph::readGraphFile(std::string filename) {
        std::cout << "\tReading file " + filename << std::endl;

        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::string line;
            int ncells, id, nvels;
            int idx = 0;
            double x, y ,z;
            std::vector<cgal_point_3> pp;
            std::vector<int> dd;
            while (getline(datafile, line)){
                if (line.size() > 1){
                    std::istringstream inp(line.c_str());
                    inp >> x;
                    inp >> y;
                    inp >> z;
                    pp.push_back(cgal_point_3(x, y, z));
                    dd.push_back(idx);
                    idx++;
                    inp >> ncells;
                    inp >> nvels;
                    std::vector<int> cells;
                    for (int i = 0; i < ncells; ++i){
                        inp >> id;
                        cells.push_back(id);
                    }
                    std::vector<int> vels;
                    for (int i = 0; i < nvels; ++i){
                        inp >> id;
                        vels.push_back(id);
                    }
                    neighElem.push_back(cells);
                    velids.push_back(vels);
                }
            }

            std::cout << neighElem.size() << ", " << velids.size() << std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            Tree.insert(boost::make_zip_iterator(boost::make_tuple( pp.begin(), dd.begin() )),
                        boost::make_zip_iterator(boost::make_tuple( pp.end(), dd.end() ) )  );
            Tree.build();
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "\tPoint Set Building time: " << elapsed.count() << std::endl;
            return true;
        }
    }

    bool CellGraph::getNearVelocities(vec3 &p, std::vector<int>& outvel) {
        cgal_point_3 query(p.x, p.y, p.z);
        K_neighbor_search_int search(Tree, query, 1);
        int idx = -9;
        for (K_neighbor_search_int::iterator it = search.begin(); it != search.end(); it++){
            idx = boost::get<1>(it->first);
        }
        if (idx == -9){
            std::cout << "idx not found" << std::endl;
            return false;
        }
        if (idx >= 0 && idx < velids.size()){
            // Add the velocities of the same cell
            for (int i = 0; i < velids[idx].size(); ++i){
                outvel.push_back(velids[idx][i]);
            }
            // Loop through the neighbor elements and add the velocity ids of those elements
            for (int i = 0; i < neighElem[idx].size(); ++i){
                int elid = neighElem[idx][i];
                if (elid >= 0 && elid < velids.size()){
                    for (int j = 0; j < velids[elid].size(); ++j){
                        outvel.push_back(velids[elid][j]);
                    }
                }
                else{
                    std::cout << "elid " << elid << " out of range " << velids.size() << std::endl;
                }
            }
            return true;
        }
        else{
            std::cout << "idx " << idx << " out of range " << velids.size() << std::endl;
            return false;
        }
    }
}


