#pragma once
#include <boost/lexical_cast.hpp>

#include "ichnos_structures.h"
#include "ichnos_utils.h"

namespace ICHNOS {

	class DomainBase {
	public:
		DomainBase(DomainOptions& Dopt_in);
		virtual void bisInPolygon(vec3& p, bool& tf) {};
		virtual void getTopBottomElevation(vec3 p, double& top, double& bottom) {};
	protected:
		multiPoly domainPoly;
		DomainOptions Dopt;
	};

	DomainBase::DomainBase(DomainOptions& Dopt_in)
		:
		Dopt(Dopt_in)
	{}


	class Domain2D : public DomainBase {
	public:
		Domain2D(DomainOptions& Dopt_in);

		void bisInPolygon(vec3& p, bool& tf);
		void getTopBottomElevation(vec3 p, double& top, double& bottom);

	private:
		cgal_Delaunay_2	Top_tria;
		cgal_Delaunay_2 Bot_tria;
		coord_map top_values;
		coord_map bot_values;
		
		float top_value;
		float bot_value;
		bool isTopConst;
		bool isBottomConst;
	};

	Domain2D::Domain2D(DomainOptions& Dopt_in)
		:
		DomainBase(Dopt_in)
	{
		domainPoly.readfromFile(Dopt.polygonFile);
		//READ::readPolygonDomain(Dopt.OutlineFile, outlinePoly);
		
		if (is_input_scalar(Dopt.TopElevationFile)) {
			isTopConst = true;
			top_value = std::stof(Dopt.TopElevationFile);
		}
		else {
			READ::read2DscalarField(Dopt.TopElevationFile, Top_tria, top_values);
			isTopConst = false;
		}

		if (is_input_scalar(Dopt.BottomeElevationFile)) {
			isBottomConst = true;
			bot_value = std::stof(Dopt.BottomeElevationFile);
		}
		else {
			READ::read2DscalarField(Dopt.BottomeElevationFile, Bot_tria, bot_values);
			isBottomConst = false;
		}
	}

	void Domain2D::bisInPolygon(vec3& p, bool& tf) {
		tf = domainPoly.is_point_in(p.x, p.y);
		//boostPoint pnt(p.x, p.y);
		//if (boost::geometry::within(pnt, outlinePoly))
		//	tf = true;
		//else
		//	tf = false;
	}

	void Domain2D::getTopBottomElevation(vec3 p, double& top, double& bottom) {
		if (!isTopConst) {
			top = interpolateScatter2D(Top_tria, top_values, p.x, p.y);
		}
		else
			top = top_value;

		if (!isBottomConst) {
			bottom = interpolateScatter2D(Bot_tria, bot_values, p.x, p.y);
		}
		else
			bottom = bot_value;
	}

}