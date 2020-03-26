#pragma once
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include <iostream>
#include <vector>

#include <nanoflann.hpp>

namespace ICHNOS {

	struct WellOptions {
		int Nparticles;
		int Nlayer;
		double Radius;
	};

	struct vec3
	{
		double x = 0;
		double y = 0;
		double z = 0;

		vec3(double x = 0, double y = 0, double z = 0)
			:
			x(x), y(y), z(z)
		{}

		vec3 operator*(const vec3& a) const {
			return vec3(a.x * x, a.y * y, a.z * z);
		}

		vec3 operator*(const double& a) const {
			return vec3(a * x, a * y, a * z);
		}

		vec3 operator+(const vec3& a) const {
			return vec3(a.x + x, a.y + y, a.z + z);
		}

		vec3 operator-(const vec3& a) const {
			//return vec3(a.x - x, a.y - y, a.z - z);
			return vec3(x - a.x, y - a.y, z - a.z);
		}

		double len() {
			return std::sqrt(x * x + y * y + z * z);
		}

		bool isZero() {
			bool tf = len() < 0.0000000001;
			return tf;
		}
		bool isInvalid() {
			bool tf = std::sqrt((x + 99999) * (x + 99999) + (y + 99999) * (y + 99999) + (z + 99999) * (z + 99999)) < 0.0000000001;
			return tf;
		}
	};

	template <typename S>
	struct pointCloud {
		std::vector<vec3> pts;
		std::vector<S> vel;
		std::vector<int> proc;

		double Radious;
		double Power;
		double Threshold;
		bool InterpolateOutside;


		// methods required by nanoflann ------
		// Must return the number of data points
		inline size_t kdtree_get_point_count() const { return pts.size(); }

		// Returns the dim'th component of the idx'th point in the class:
		// Since this is inlined and the "dim" argument is typically an immediate value, the
		//  "if/else's" are actually solved at compile time.
		inline double kdtree_get_pt(const size_t idx, const size_t dim) const
		{
			if (dim == 0) return pts[idx].x;
			else if (dim == 1) return pts[idx].y;
			else return pts[idx].z;
		}

		inline vec3 kdtree_get_pnt_vec(const size_t idx) const {
			return pts[idx];
		}

		inline S kdtree_get_vel_vec(const size_t idx) const {
			return vel[idx];
		}
		inline int kdtree_get_proc(const size_t idx) const {
			return proc[idx];
		}

		// Optional bounding-box computation: return false to default to a standard bbox computation loop.
		//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
		//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
		template <class BBOX>
		bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
	};

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
		FAR_AWAY
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
		it = ExitReasonsMap.find(er);
		if (it != ExitReasonsMap.end())
			return it->second;
		else {
			return "INVALID";
		}
	}

	struct ParticleOptions {
		SolutionMethods method;
		int StuckIterations;
		double MaxStepSize;
		double MinStepSize;
		double StepSize;
		double ToleranceStepSize;
		double Direction;
		int MaxIterationsPerStreamline;
		int MaxProcessorExchanges;
		std::string ParticleFile;
		std::string WellFile;
		int ParticlesInParallel;
		std::string OutputFile;
	};

	struct DomainOptions {
		std::string OutlineFile;
		std::string TopElevationFile;
		std::string BottomeElevationFile;
	};


	class Particle {
	public:
		Particle(vec3 p);
		Particle(vec3 p, vec3 v, int pid, int proc);
		Particle(vec3 p, vec3 v, int proc,  Particle prev);
		vec3 getP() const { return p; }
		vec3 getV() const { return v; }
		double getAge() const { return age; }
		void SetAge(double a) { age = a; }
		void SetV(vec3 vel) { v = vel; }
		int getPid() const { return pid; }
		void displayAsVEX(bool printAttrib);
		void SetProc(int procId) { proc = procId; }
		int getProc() { return proc; }
	private:
		int pid;
		vec3 p;
		vec3 v;
		double age;
		int proc;
	};

	Particle::Particle(vec3 p)
		:
		p(p)
	{
		pid = 0;
		age = 0;
		proc = -9;
	}

	Particle::Particle(vec3 p, vec3 v, int pid, int proc)
		:
		pid(pid),
		p(p),
		v(v),
		proc(proc)
	{
		age = 0;
	}

	Particle::Particle(vec3 p, vec3 v, int proc, Particle prev)
		:
		p(p),
		v(v),
		proc(proc)
	{
		pid = prev.getPid() + 1;
		double d = (p - prev.getP()).len();
		age = prev.getAge() + d / prev.getV().len();
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
		Streamline(int Eid_in, int Sid_in, Particle particle, vec3 BL, vec3 BU, int N);
		void InitialVelocity(vec3 v) {
			SL.back().SetV(v);
		}
		void AddParticle(Particle p);
		Particle getLastParticle() { return SL.back(); }
		void Close(ExitReason ER);
		int StuckIter() { return IterationsWithoutExpansion; }
		int getEid() { return Eid; }
		int getSid() { return Sid; }
		vec3 getBBlow() { return BL; }
		vec3 getBBupp() { return BU; }
		unsigned int size() { return SL.size(); }
		Particle getParticle(int i) { return SL.at(i); }


	private:
		int Eid;
		int Sid;
		std::vector<Particle > SL;
		ExitReason exitreason;
		vec3 BL;
		vec3 BU;
		int IterationsWithoutExpansion = 0;
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
	}

	Streamline::Streamline(int Eid_in, int Sid_in, Particle particle, vec3 BL, vec3 BU, int N)
		:
		Eid(Eid_in),
		Sid(Sid_in),
		BL(BL),
		BU(BU),
		IterationsWithoutExpansion(N)
	{
		SL.clear();
		SL.push_back(particle);
	}

	void Streamline::AddParticle(Particle p) {
		SL.push_back(p);
		vec3 pp = p.getP();
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
		RKcoef.push_back(std::vector<double>{1 / 4, 1 / 4});
		RKcoef.push_back(std::vector<double>{3 / 8, 3 / 32, 9 / 32});
		RKcoef.push_back(std::vector<double>{12 / 13, 1932 / 2197, -7200 / 2197, 7296 / 2197});
		RKcoef.push_back(std::vector<double>{1, 439 / 216, -8, 3680 / 513, -845 / 4104});
		RKcoef.push_back(std::vector<double>{1/2, -8 / 27, 2, -3544 / 2565, 1859 / 4104, -11 / 40});
		RKcoef.push_back(std::vector<double>{1, 25 / 216, 1408/2565, 2197 / 4101, -1 / 5});
		RKcoef.push_back(std::vector<double>{1, 16 / 135, 6656 / 12825, 28561 / 56430, -9 / 50, 2/55});

		return RKcoef;
	}
}

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, ICHNOS::pointCloud< ICHNOS::vec3 > >, ICHNOS::pointCloud<ICHNOS::vec3 >, 3 > nano_kd_tree_vector;
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, ICHNOS::pointCloud< double > >, ICHNOS::pointCloud<double >, 3 > nano_kd_tree_scalar;

typedef boost::geometry::model::d2::point_xy<double> boostPoint;
typedef boost::geometry::model::polygon<boostPoint> boostPolygon;

