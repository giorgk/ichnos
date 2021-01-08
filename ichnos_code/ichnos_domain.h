#pragma once
#include <boost/lexical_cast.hpp>

#include "ichnos_structures.h"
#include "ichnos_utils.h"

namespace ICHNOS {

	class DomainBase {
	public:
		DomainBase(DomainOptions& Dopt_in);
		virtual void bisInPolygon(vec3& p, bool& tf) {};
		virtual void bisInProcessorPolygon(vec3 p, bool& tf) {};
		virtual void bisInExpandedPolygon(vec3 p, bool& tf) {};
		virtual void getTopBottomElevation(vec3 p, double& top, double& bottom) {};
		boostPolygon getProcessorDomain();
	protected:
		multiPoly domainPoly;
		DomainOptions Dopt;
		boostPolygon ProcessorDomain;
		boostPolygon ExpandedDomain;
	};

	DomainBase::DomainBase(DomainOptions& Dopt_in)
		:
		Dopt(Dopt_in)
	{}

	boostPolygon DomainBase::getProcessorDomain() {
		return ProcessorDomain;
	}

	class Domain2D : public DomainBase {
	public:
		Domain2D(DomainOptions& Dopt_in);

		void bisInPolygon(vec3& p, bool& tf);
		void bisInProcessorPolygon(vec3 p, bool& tf);
		void bisInExpandedPolygon(vec3 p, bool& tf);
		void getTopBottomElevation(vec3 p, double& top, double& bottom);

	private:
		PointSet2 TopSet;
		PointSet2 BotSet;
		//cgal_Delaunay_2	Top_tria;
		//cgal_Delaunay_2 Bot_tria;
		//coord_map top_values;
		//coord_map bot_values;
		
		double top_value;
		double bot_value;
		bool isTopConst;
		bool isBottomConst;
		bool useSamePoints = false;
		double TopRadius = 1000;
		double BotRadius = 1000;
		double TopPower = 3;
		double BotPower = 3;

		void interpolateSets(vec3& p, double& top, double& bot);
	};

	Domain2D::Domain2D(DomainOptions& Dopt_in)
		:
		DomainBase(Dopt_in)
	{
		domainPoly.readfromFile(Dopt.polygonFile);
		
		READ::readProcessorDomain(Dopt.processorDomainFile, ProcessorDomain, /*dbg_rank*/ Dopt.myRank);
		READ::readProcessorDomain(Dopt.expandedDomainFile, ExpandedDomain, Dopt.myRank);

		
		//READ::readPolygonDomain(Dopt.OutlineFile, outlinePoly);
		bool getTopFromFile = false;
		bool getBotFromFile = false;
		
		if (is_input_scalar(Dopt.TopElevationFile)) {
			isTopConst = true;
			top_value = std::stof(Dopt.TopElevationFile);
			getTopFromFile = false;
			useSamePoints = false;
		}
		else {
			getTopFromFile = true;
			TopRadius = Dopt.TopRadius;
			TopPower = Dopt.TopPower;
			top_value = 0;
			isTopConst = false;
		}

		if (is_input_scalar(Dopt.BottomeElevationFile)) {
			isBottomConst = true;
			bot_value = std::stof(Dopt.BottomeElevationFile);
			getBotFromFile = false;
		}
		else{
			if (Dopt.BottomeElevationFile.empty()){
				if (is_input_scalar(Dopt.TopElevationFile)){
					std::cout << "You Can't have constant Top elevation and Empty bottom elevation" << std::endl;
				}
				else{
					getBotFromFile = false;
					useSamePoints = true;
				}
			}
			else{
				getBotFromFile = true;
				BotPower = Dopt.BotPower;
				BotRadius = Dopt.BotRadius;
			}
			bot_value = 0;
			isBottomConst = false;
		}

		if (getTopFromFile){
			if (useSamePoints)
				READ::readTopBot(Dopt.TopElevationFile, TopSet, true, true);
			else
				READ::readTopBot(Dopt.TopElevationFile, TopSet, true, false);
		}
		if (getBotFromFile){
			READ::readTopBot(Dopt.BottomeElevationFile, BotSet,  false, true);
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

	void Domain2D::bisInProcessorPolygon(vec3 p, bool& tf) {
		tf = boost::geometry::within(boostPoint(p.x, p.y), ProcessorDomain);
	}

	void Domain2D::bisInExpandedPolygon(vec3 p, bool& tf) {
		tf = boost::geometry::within(boostPoint(p.x, p.y), ExpandedDomain);
	}

	void Domain2D::getTopBottomElevation(vec3 p, double& top, double& bottom) {
		interpolateSets(p, top, bottom);
	}

	void Domain2D::interpolateSets(vec3& p, double& top, double& bot){
		std::list<Vertex_handle2D> LV;
		cgal_point_2 cntr(p.x, p.y);
		double dist, w, sumW, sumWT, sumWB, xx, yy;
		if (useSamePoints){
			CGAL::Circle_2<K> rc(cntr, TopRadius);
			TopSet.range_search(rc,std::back_inserter(LV));
			std::list<Vertex_handle2D>::const_iterator it = LV.begin();
			sumWT = 0;
			sumWB = 0;
			sumW = 0;
			for (;it != LV.end(); ++it){
				xx = (*it)->point().x() - p.x;
				yy = (*it)->point().y() - p.y;
				dist = sqrt(xx*xx + yy*yy);
				if (dist < 0.01){
					top = (*it)->info().top;
					bot = (*it)->info().bot;
					return;
				}
				w = 1/std::pow(dist,TopPower);
				sumW += w;
				sumWT += w*(*it)->info().top;
				sumWB += w*(*it)->info().bot;
			}
			top = sumWT/sumW;
			bot = sumWB/sumW;
		}
		else{
			if (!isTopConst){
				LV.clear();
				CGAL::Circle_2<K> rc_top(cntr, TopRadius);
				TopSet.range_search(rc_top,std::back_inserter(LV));
				std::list<Vertex_handle2D>::const_iterator it = LV.begin();
				sumWT = 0;
				sumW = 0;
				bool calc_mean = true;
				for (;it != LV.end(); ++it){
					xx = (*it)->point().x() - p.x;
					yy = (*it)->point().y() - p.y;
					dist = sqrt(xx*xx + yy*yy);
					if (dist < 0.01){
						top = (*it)->info().top;
						calc_mean = false;
						break;
					}
					w = 1/std::pow(dist,TopPower);
					sumW += w;
					sumWT += w*(*it)->info().top;
				}
				if (calc_mean)
					top = sumWT/sumW;

			}
			else{
				top = top_value;
			}
			if (!isBottomConst){
				LV.clear();
				CGAL::Circle_2<K> rc_bot(cntr, BotRadius);
				BotSet.range_search(rc_bot, std::back_inserter(LV));
				std::list<Vertex_handle2D>::const_iterator it = LV.begin();
				sumWB = 0;
				sumW = 0;
				bool calc_mean = true;
				for (;it != LV.end(); ++it){
					xx = (*it)->point().x() - p.x;
					yy = (*it)->point().y() - p.y;
					dist = sqrt(xx*xx + yy*yy);
					if (dist < 0.01){
						bot = (*it)->info().bot;
						break;
					}
					w = 1/std::pow(dist,TopPower);
					sumW += w;
					sumWB += w*(*it)->info().bot;
				}
				if (calc_mean)
					bot = sumWB/sumW;
			}
			else{
				bot = bot_value;
			}
		}
	}
}