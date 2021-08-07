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


	struct NPSAT_data {
		int proc = -9;
		int id = -9;
		double diameter = 0;
		double ratio = 0;
		vec3 v;
	};

	struct elev_data{
		double top = 0;
		double bot = 0;
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
	typedef boost::tuple<cgal_point_3, NPSAT_data> pnt_int;
	typedef CGAL::Search_traits_3<K> Traits_base;
	typedef CGAL::Search_traits_adapter<pnt_int, CGAL::Nth_of_tuple_property_map<0,pnt_int>, Traits_base> search_traits;
	typedef CGAL::Fuzzy_iso_box<search_traits> Fuzzy_iso_box;
	typedef CGAL::Kd_tree<search_traits> search_tree;

	

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


    enum class coordDim {vx, vy, vz};
    enum class TimeInterpType {NEAREST, LINEAR};

    struct TRANS_data{
        int proc = -9;
        int id = -9;
        double diameter = 0;
        double ratio = 0;
    };

    // Stochastic tree
    typedef boost::tuple<cgal_point_3, TRANS_data> pnt_trans;
    typedef CGAL::Search_traits_adapter<pnt_trans, CGAL::Nth_of_tuple_property_map<0, pnt_trans>, Traits_base> search_traits_trans;
    typedef CGAL::Fuzzy_iso_box<search_traits_trans> Fuzzy_iso_box_trans;
    typedef CGAL::Kd_tree<search_traits_trans> search_tree_trans;

    class VelTR{
    public:
        VelTR();
        void init(int np, int nt);
        void setVELvalue(double v, int pnt, int step, coordDim dim);
        void setTSvalue(double v, int step);
        void setTSvalue(std::vector<double>& TS_in);
        void findIIT(double x, int &i1, int &i2, double &t);
        vec3 getVelocity(int pnt, int i1, int i2, double t, TimeInterpType tp);
    private:
        std::vector<std::vector<double>> VX;
        std::vector<std::vector<double>> VY;
        std::vector<std::vector<double>> VZ;

        std::vector<double> TS;

        int nPoints;
        int nSteps;
        void findTimeStepIndex(int &i, int &ii, double x);

    };

    VelTR::VelTR()
        :
        nPoints(-9),
        nSteps(-9)
    {}

    void VelTR::init(int np, int nt) {
        nPoints = np;
        nSteps = nt;
        VX.clear();
        VY.clear();
        VZ.clear();
        VX.resize(nPoints, std::vector<double>(nSteps, 0));
        VY.resize(nPoints, std::vector<double>(nSteps, 0));
        VZ.resize(nPoints, std::vector<double>(nSteps, 0));
        TS.clear();
        TS.resize(nSteps, 0);
    }

    void VelTR::setVELvalue(double v, int pnt, int step, coordDim dim){
        if (pnt >= nPoints || step >= nSteps || pnt < 0  || step < 0)
            return;

        if (dim == coordDim::vx){
            VX[pnt][step] = v;
        }
        else if (dim == coordDim::vy){
            VY[pnt][step] = v;
        }
        else if (dim == coordDim::vz){
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

    void VelTR::findIIT(double x, int &i1, int &i2, double &t){
        if (x <= TS[0]){
            i1 = 0;
            i2 = 0;
            t = 0.0;
            return;
        }
        if (x >= TS[TS.size()-1]){
            i1 = TS.size()-1;
            i2 = TS.size()-1;
            t = 1.0;
            return;
        }
        i1 = 0;
        i2 = TS.size()-1;
        findTimeStepIndex(i1, i2, x);
        t = (x - TS[i1]) / (TS[i2] - TS[i1]);
    }

    vec3 VelTR::getVelocity(int pnt, int i1, int i2, double t, TimeInterpType tp){
        if (tp == TimeInterpType::LINEAR){
            vec3 v1(VX[pnt][i1], VY[pnt][i1], VZ[pnt][i1]);
            vec3 v2(VX[pnt][i2], VY[pnt][i2], VZ[pnt][i2]);
            return v1*(1-t) + v2*t;

        }
        else{
            return vec3(VX[pnt][i1], VY[pnt][i1], VZ[pnt][i1]);
        }
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
		it = ExitReasonMap.find(exitreason);
		if (it != ExitReasonMap.end())
			return it->second;
		else {
			return ExitReason::NO_REASON;
		}
	}


	struct ParticleOptions {
		
		SolutionMethods method;
		// The maximum number of iteratiosn that each streamline is alloweded not to expand
		int StuckIterations;
		// When the method is adaptive this is the maximum step size
		double MaxStepSize;
		// When the method is adaptive this is the minimum step size
		double MinStepSize;
		// This is the step size of the initial step size for the adaptive
		double StepSize;
		// This is the minimum step size at the exit side of the particles
		double minExitStepSize;
		//When the method is adaptive this is the tolerance
		double ToleranceStepSize;
		// When the method is adaptive this modifies the rate of change of the step size.
		double increasRatechange;

		// When the method is adaptive this limits the upper threshold of step decrease.
		// For example if the decrease is lower that 1 but very close to one then it may be possible that
		// the algorithm will repeat the decreasing of the step size many times by a very small amount.
		double limitUpperDecreaseStep;

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
		// During gathering puut all streamlines into one file.
		// If false each file will contain ParticlesInParallel streamlines per file
		bool gatherOneFile;
		std::string configfile;
		// Number of realizations
		int Nrealizations;

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
		std::string expandedDomainFile;
		int myRank;
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
		std::cout << "p = addpoint(0,{" << p.x << "," << p.z << "," << p.y << ",});";
		if (printAttrib)
			std::cout << " setpointattrib(0,'N',p,{" << vn.x << "," << vn.z << "," << vn.y << "},'set');";
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


	private:
		int Eid;
		int Sid;
		std::vector<Particle > SL;
		ExitReason exitreason = ExitReason::NO_REASON;
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
		//double age;
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
		}

		boost::random::uniform_real_distribution<float> uniformDistribution;
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
}


