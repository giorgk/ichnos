#pragma once
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <iostream>
#include <fstream>
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

		vec3 normalize() {
			double length = len();
			return (vec3(x / length, y / length, z / length));
		}
		void zero() {
			x = 0.0;
			y = 0.0;
			z = 0.0;
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
		int NmaxPnts;
		int NminPnts;


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

	enum class interpType { INGORE, SCALAR, CLOUD };

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
		// THe maximum number of iteratiosn that each streamline is alloweded not to expand
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

		// The direction of the particle tracing
		double Direction;
		// The maximum steps for each streamline
		int MaxIterationsPerStreamline;
		// The maximum number of exchanges between processors
		int MaxProcessorExchanges;
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
	};

	struct DomainOptions {
		std::string polygonFile;
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
		//double getAge() const { return age; }
		//void SetAge(double a) { age = a; }
		void SetV(vec3 vel) { v = vel; }
		int getPid() const { return pid; }
		void displayAsVEX(bool printAttrib);
		void SetProc(int procId) { proc = procId; }
		int getProc() { return proc; }
	private:
		int pid;
		vec3 p;
		vec3 v;
		//double age;
		int proc;
	};

	Particle::Particle(vec3 p)
		:
		p(p)
	{
		pid = 0;
		//age = 0;
		proc = -9;
	}

	Particle::Particle(vec3 p, vec3 v, int pid, int proc)
		:
		pid(pid),
		p(p),
		v(v),
		proc(proc)
	{
		//age = 0;
	}

	Particle::Particle(vec3 p, vec3 v, int proc, Particle prev)
		:
		p(p),
		v(v),
		proc(proc)
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
		ExitReason exitreason = ExitReason::NO_REASON;
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

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, ICHNOS::pointCloud< ICHNOS::vec3 > >, ICHNOS::pointCloud<ICHNOS::vec3 >, 3 > nano_kd_tree_vector;
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, ICHNOS::pointCloud< double > >, ICHNOS::pointCloud<double >, 3 > nano_kd_tree_scalar;
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, ICHNOS::pointCloud< int > >, ICHNOS::pointCloud<int >, 3 > nano_kd_tree_int;

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

}


