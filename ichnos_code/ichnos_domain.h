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
		boostPolygon outlinePoly;
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
		pointCloud<double> TopElevation;
		pointCloud<double> BottomElevation;
		std::unique_ptr<nano_kd_tree_scalar> TopTree;
		std::unique_ptr<nano_kd_tree_scalar> BottomTree;
		
		float top_value;
		float bot_value;
		bool isTopCloud;
		bool isBottomCloud;
		
	};

	Domain2D::Domain2D(DomainOptions& Dopt_in)
		:
		DomainBase(Dopt_in)
	{
		READ::readPolygonDomain(Dopt.OutlineFile, outlinePoly);
		
		if (is_input_scalar(Dopt.TopElevationFile)) {
			isTopCloud = false;
			top_value = std::stof(Dopt.TopElevationFile);
		}
		else {
			READ::read2DscalarField(Dopt.TopElevationFile, TopElevation);
			this->TopTree = std::unique_ptr<nano_kd_tree_scalar>(new nano_kd_tree_scalar(3, TopElevation, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
			TopTree->buildIndex();
			isTopCloud = true;
		}

		if (is_input_scalar(Dopt.BottomeElevationFile)) {
			isBottomCloud = false;
			bot_value = std::stof(Dopt.BottomeElevationFile);
		}
		else {
			READ::read2DscalarField(Dopt.BottomeElevationFile, BottomElevation);
			this->BottomTree = std::unique_ptr<nano_kd_tree_scalar>(new nano_kd_tree_scalar(3, BottomElevation, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
			BottomTree->buildIndex();
			isBottomCloud = true;
		}
	}

	void Domain2D::bisInPolygon(vec3& p, bool& tf) {
		boostPoint pnt(p.x, p.y);
		if (boost::geometry::within(pnt, outlinePoly))
			tf = true;
		else
			tf = false;
	}

	void Domain2D::getTopBottomElevation(vec3 p, double& top, double& bottom) {
		if (isTopCloud)
			top = interpolateScalarTree(TopTree, p);
		else
			top = top_value;

		if (isBottomCloud)
			bottom = interpolateScalarTree(BottomTree, p);
		else
			bottom = bot_value;
	}

}