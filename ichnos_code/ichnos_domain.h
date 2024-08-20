#pragma once
#include <boost/lexical_cast.hpp>

//#include "gridInterp.h"

#include "ichnos_structures.h"
#include "ichnos_utils.h"
#include "ichnos_mesh.h"

namespace ICHNOS {

	class DomainBase {
	public:
		DomainBase(DomainOptions& Dopt_in);
		virtual void bisInPolygon(vec3& p, bool& tf) {};
		virtual void bisInProcessorPolygon(vec3 p, bool& tf) {};
		//virtual void bisInExpandedPolygon(vec3 p, bool& tf) {};
		virtual void bisInNearAttractor(vec3 p, bool& tf) {};
		virtual void getTopBottomElevation(vec3 p, double& top, double& bottom, bool& tf) {};
		boostPolygon getProcessorDomain();
	protected:
		multiPoly domainPoly;
		DomainOptions Dopt;
		boostPolygon ProcessorDomain;
		//boostPolygon ExpandedDomain;
		int nProc;
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

        enum type{
            CONST,
            CLOUD,
            MESH2D,
            GRID
        };

		void bisInPolygon(vec3& p, bool& tf);
		void bisInProcessorPolygon(vec3 p, bool& tf);
		//void bisInExpandedPolygon(vec3 p, bool& tf);
		void bisInNearAttractor(vec3 p, bool& tf);
		void getTopBottomElevation(vec3 p, double& top, double& bottom, bool& tf);

	private:
		PointSet2 TopSet;
		PointSet2 BotSet;
		PointSet2 AttrctSet;
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
		bool useAttractors = false;
		double AttractRadius = 10;
		Mesh2D TopMesh;
        Mesh2D BotMesh;
        std::vector<double> TopNodes;
        std::vector<double> BotNodes;
        type Top_interpolation_type;
        type Bot_interpolation_type;

        double invTransfTol = 0.001;

        //GRID_INTERP::interp<2> gridTop;
        //GRID_INTERP::interp<2> gridBot;


		bool interpolateSets(vec3& p, double& top, double& bot);
		bool readTopBot(std::string filename, bool top, bool bot);
		bool CloudInterpolate(vec3& p, double& top, double& bot, bool btop, bool bbot);
		bool MeshInterpolate(vec3& p, double& top, double& bot, bool btop, bool bbot);
	};

	Domain2D::Domain2D(DomainOptions& Dopt_in)
		:
		DomainBase(Dopt_in)
	{
		domainPoly.readfromFile(Dopt.polygonFile);
		nProc = Dopt_in.nProc;

		if (nProc > 1 & !Dopt_in.RunAsThread){
            //std::cout << "Dopt.myRank " << Dopt.myRank << std::endl;
            READ::readProcessorDomain(Dopt.processorDomainFile, ProcessorDomain, /*dbg_rank*/ Dopt.myRank);
            //READ::readProcessorDomain(Dopt.expandedDomainFile, ExpandedDomain, Dopt.myRank);
		}


		
		//READ::readPolygonDomain(Dopt.OutlineFile, outlinePoly);
		bool getTopFromFile = false;
		bool getBotFromFile = false;
		
		if (is_input_scalar(Dopt.TopElevationFile)) {
			Top_interpolation_type = type::CONST;
			top_value = std::stof(Dopt.TopElevationFile);
			getTopFromFile = false;
			useSamePoints = false;
		}
		else {
			getTopFromFile = true;
			top_value = 0;
		}

		if (is_input_scalar(Dopt.BottomeElevationFile)) {
			Bot_interpolation_type = type::CONST;
			bot_value = std::stof(Dopt.BottomeElevationFile);
			getBotFromFile = false;
		}
		else{
			if (Dopt.BottomeElevationFile.empty()){
				if (is_input_scalar(Dopt.TopElevationFile)){
					std::cout << "You Can't have constant Top elevation and empty bottom elevation" << std::endl;
				}
				else{
					getBotFromFile = false;
					useSamePoints = true;
				}
			}
			else{
				getBotFromFile = true;
			}
			bot_value = 0;
		}

		if (getTopFromFile){
			if (useSamePoints) {
			    readTopBot(Dopt.TopElevationFile, true, true);
			}
			else{
                readTopBot(Dopt.TopElevationFile, true, false);
			}
		}
		if (getBotFromFile){
            readTopBot(Dopt.BottomeElevationFile, false, true);
		}

		if (!Dopt.AttractorsFile.empty()) {
			useAttractors = true;
			READ::readTopBot(Dopt.AttractorsFile, AttrctSet, true, true);
			AttractRadius = Dopt.AttractRadius;
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
        if (nProc == 1 || Dopt.RunAsThread){
            tf = true;
            return;
        }
		tf = boost::geometry::within(boostPoint(p.x, p.y), ProcessorDomain);
	}

	/*
	void Domain2D::bisInExpandedPolygon(vec3 p, bool& tf) {
        if (nProc == 1){
            tf = true;
            return;
        }
		tf = boost::geometry::within(boostPoint(p.x, p.y), ExpandedDomain);
	}
	*/

	void Domain2D::getTopBottomElevation(vec3 p, double& top, double& bottom, bool& tf) {
		tf = interpolateSets(p, top, bottom);
	}

	void Domain2D::bisInNearAttractor(vec3 p, bool& tf) {
		tf = false;
		if (!useAttractors) {
			return;
		}
		std::list<Vertex_handle2D> LV;
		cgal_point_2 cntr(p.x, p.y);
		CGAL::Circle_2<K> rc(cntr, AttractRadius);
		AttrctSet.range_search(rc, std::back_inserter(LV));
		if (!LV.empty()) {
			std::list<Vertex_handle2D>::const_iterator it = LV.begin();
			for (; it != LV.end(); ++it) {
				double bot_attract = (*it)->info().bot - AttractRadius;
				double top_attract = (*it)->info().top + AttractRadius;
				if (p.z > bot_attract && p.z < top_attract) {
					tf = true;
					break;
				}
			}
		}
	}

	bool Domain2D::interpolateSets(vec3& p, double& top, double& bot){
	    if (useSamePoints){
	        if (Top_interpolation_type == type::CLOUD){
                return CloudInterpolate(p, top, bot, true, true);
	        }
	        else if (Top_interpolation_type == type::MESH2D){
                return MeshInterpolate(p, top, bot, true, true);
	        }
	    }
	    else{
	        bool top_tf = false;
	        bool bot_tf = false;
	        double dummy;
            switch (Top_interpolation_type) {
                case type::CONST:
                {
                    top = top_value;
                    top_tf = true;
                    break;
                }
                case type::CLOUD:
                {
                    top_tf = CloudInterpolate(p, top, dummy, true, false);
                    break;
                }
                case type::MESH2D:
                {
                    top_tf = MeshInterpolate(p, top, dummy, true, false);
                    break;
                }
            }
            switch (Bot_interpolation_type) {
                case type::CONST:
                {
                    bot = bot_value;
                    bot_tf = true;
                    break;
                }
                case type::CLOUD:
                {
                    bot_tf = CloudInterpolate(p, dummy, bot, false, true);
                    break;
                }
                case type::MESH2D:
                {
                    bot_tf = MeshInterpolate(p, dummy, bot, false, true);
                    break;
                }
            }
            if (top_tf && bot_tf)
                return true;
            return false;
	    }
	    return false;
	}

	bool Domain2D::readTopBot(std::string filename, bool top, bool bot) {
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file" << filename << std::endl;
            return false;
        }
        else{
            std::vector< std::pair<cgal_point_2,elev_data> > xy_data;
            type interpolation_type;
            std::string line;
            { // Read the TYPE
                getline(datafile, line);
                if (line.compare("CLOUD") == 0){
                    interpolation_type = type::CLOUD;
                }
                else if (line.compare("MESH2D") == 0){
                    interpolation_type = type::MESH2D;
                }
                else if (line.compare("GRID") == 0){
                    interpolation_type = type::GRID;
                }
                else{
                    std::cout << "[" << line << "] is unknown interpolation type for the top/bottom function" << std::endl;
                    return false;
                }
            }

            if (top)
                Top_interpolation_type = interpolation_type;
            if (bot)
                Bot_interpolation_type = interpolation_type;

            {// read data
                getline(datafile, line);
                std::istringstream inp(line.c_str());
                if (interpolation_type == type::CLOUD){
                    double rad, pw;
                    inp >> rad;
                    inp >> pw;
                    if (top){
                        TopRadius = rad*rad;
                        TopPower = pw;
                    }
                    if (bot){
                        BotRadius = rad*rad;
                        BotPower = pw;
                    }
                    double x, y;
                    while (getline(datafile, line)){
                        if (line.size() > 1){
                            elev_data el_data;
                            std::istringstream inp1(line.c_str());
                            inp1 >> x;
                            inp1 >> y;
                            if (top && bot){
                                inp1 >> el_data.top;
                                inp1 >> el_data.bot;
                                if (el_data.top < el_data.bot) {
                                    std::cerr << "Point:(" << x << "," << y << ") has lower top" << std::endl;
                                }
                            }
                            else if (top){
                                inp1 >> el_data.top;
                            }
                            else if (bot){
                                inp1 >> el_data.bot;
                            }
                            cgal_point_2 p(x, y);
                            xy_data.push_back(std::make_pair(p, el_data));
                        }
                    }
                    if ((!useSamePoints && top) || useSamePoints){
                        TopSet.insert(xy_data.begin(), xy_data.end());
                    }
                    else if (!useSamePoints && bot){
                        BotSet.insert(xy_data.begin(), xy_data.end());
                    }
                }
                else if (interpolation_type == type::MESH2D){
                    elev_data el_data;
                    int Nnodes, Nelem;
                    int splitQuads = 0;
                    double rad;
                    inp >> Nnodes;
                    inp >> Nelem;
                    //inp >> splitQuads;
                    inp >> rad;
                    if (top){
                        TopRadius = rad*rad;
                    }
                    if (bot){
                        BotRadius = rad*rad;
                    }
                    double x,y,v;
                    std::vector<vec3> nodes;
                    for(int i = 0; i < Nnodes; ++i){
                        getline(datafile, line);
                        std::istringstream inp1(line.c_str());
                        inp1 >> x;
                        inp1 >> y;
                        nodes.push_back(vec3(x,y,0.0));
                        if (top){
                            inp1 >> v;
                            TopNodes.push_back(v);
                        }
                        if (bot){
                            inp1 >> v;
                            BotNodes.push_back(v);
                        }
                    }
                    if ((!useSamePoints && top) || useSamePoints){
                        TopMesh.setNodes(nodes);
                        TopMesh.setInverseTransformTolerance(invTransfTol);
                    }
                    if (!useSamePoints && bot){
                        BotMesh.setNodes(nodes);
                        BotMesh.setInverseTransformTolerance(invTransfTol);
                    }
                    int idx;
                    std::vector<std::vector<int>> m;
                    for (int i = 0; i < Nelem; ++i){
                        getline(datafile, line);
                        std::istringstream inp1(line.c_str());
                        std::vector<int> tmp;
                        for (int j = 0; j < 4; ++j){
                            inp1 >> idx;
                            if (idx != 0){
                                idx = idx-1;
                                if (idx < 0){
                                    idx = idx + 1;
                                    std::cout << "I found a mesh index for top/Bottom function equal to " << idx << std::endl;
                                    return false;
                                }
                                tmp.push_back(idx);
                            }
                        }
                        if (splitQuads == 1 && tmp.size() == 4){
                            std::vector<int> tmp1{tmp[0], tmp[1], tmp[2]};
                            std::vector<int> tmp2{tmp[0], tmp[2], tmp[3]};
                            m.push_back(tmp1);
                            m.push_back(tmp2);
                        }
                        else{
                            m.push_back(tmp);
                        }
                    }
                    std::vector<vec3> cc;
                    if ((!useSamePoints && top) || useSamePoints){
                        TopMesh.setMesh(m, true);
                        TopMesh.CalculateCentroid(cc);
                        for (unsigned int i = 0;i < cc.size(); ++i){
                            cgal_point_2 p(cc[i].x, cc[i].y);
                            el_data.id = i;
                            xy_data.push_back(std::make_pair(p, el_data));
                        }
                        TopSet.insert(xy_data.begin(), xy_data.end());
                    }
                    else if (!useSamePoints && bot){
                        BotMesh.setMesh(m, true);
                        BotMesh.CalculateCentroid(cc);
                        for (unsigned int i = 0;i < cc.size(); ++i){
                            cgal_point_2 p(cc[i].x, cc[i].y);
                            el_data.id = i;
                            xy_data.push_back(std::make_pair(p, el_data));
                        }
                        BotSet.insert(xy_data.begin(), xy_data.end());
                    }
                }
                else if (interpolation_type == type::GRID){
                    std::string gridTopfFle, gridBotFile;
                    inp >> gridTopfFle;
                    inp >> gridBotFile;
                    //gridTop.getDataFromFile(gridTopfFle);
                    //gridBot.getDataFromFile(gridBotFile);
                }
            }
            datafile.close();
            return true;
        }
	}
	bool Domain2D::CloudInterpolate(vec3 &p, double &top, double &bot, bool btop, bool bbot) {
	    bool out = false;
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
            bool calc_mean = true;
            for (;it != LV.end(); ++it){
                xx = (*it)->point().x() - p.x;
                yy = (*it)->point().y() - p.y;
                dist = sqrt(xx*xx + yy*yy);
                if (dist < 0.01){
                    top = (*it)->info().top;
                    bot = (*it)->info().bot;
                    calc_mean = false;
                    out = true;
                }
                w = 1/std::pow(dist,TopPower);
                sumW += w;
                sumWT += w*(*it)->info().top;
                sumWB += w*(*it)->info().bot;
            }
            if (calc_mean) {
                top = sumWT / sumW;
                bot = sumWB / sumW;
                out = true;
            }
        }
        else{
            if (btop){
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
                        out = true;
                        break;
                    }
                    w = 1/std::pow(dist,TopPower);
                    sumW += w;
                    sumWT += w*(*it)->info().top;
                }
                if (calc_mean) {
                    top = sumWT / sumW;
                    out = true;
                }
            }

            if (bbot){
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
                        out = true;
                        break;
                    }
                    w = 1/std::pow(dist,TopPower);
                    sumW += w;
                    sumWB += w*(*it)->info().bot;
                }
                if (calc_mean){
                    bot = sumWB/sumW;
                    out = true;
                }
            }
        }
        return out;
	}

	bool Domain2D::MeshInterpolate(vec3 &p, double &top, double &bot, bool btop, bool bbot) {
	    bool out = false;
        std::list<Vertex_handle2D> LV;
        cgal_point_2 cntr(p.x, p.y);
        vec3 uv;
        if (useSamePoints){
            CGAL::Circle_2<K> rc(cntr, TopRadius);
            TopSet.range_search(rc,std::back_inserter(LV));
            std::list<Vertex_handle2D>::const_iterator it = LV.begin();
            for (;it != LV.end(); ++it){
                int idx = (*it)->info().id;
                bool tf = TopMesh.isElementInMesh(idx, p, uv);
                if (tf){
                    std::vector<int> ids;
                    TopMesh.getElementIDs(idx, ids);
                    if (ids.size() == 3){
                        top = TopNodes[ids[0]] * uv.x +
                              TopNodes[ids[1]] * uv.y +
                              TopNodes[ids[2]] * uv.z;

                        bot = BotNodes[ids[0]] * uv.x +
                              BotNodes[ids[1]] * uv.y +
                              BotNodes[ids[2]] * uv.z;
                        out = true;
                    }
                    else if (ids.size() == 4){
                        double N1, N2, N3, N4;
                        QuadShapeFunctions(uv.x, uv.y, N1, N2, N3, N4);
                        top = TopNodes[ids[0]] * N1 + TopNodes[ids[1]] * N2 +
                              TopNodes[ids[2]] * N3 + TopNodes[ids[3]] * N4;
                        bot = BotNodes[ids[0]] * N1 + BotNodes[ids[1]] * N2 +
                              BotNodes[ids[2]] * N3 + BotNodes[ids[3]] * N4;
                        out = true;
                    }
                    break;
                }
            }
        }
        else{
            if (btop){
                CGAL::Circle_2<K> rc(cntr, TopRadius);
                TopSet.range_search(rc,std::back_inserter(LV));
                std::list<Vertex_handle2D>::const_iterator it = LV.begin();
                for (;it != LV.end(); ++it){
                    int idx = (*it)->info().id;
                    bool tf = TopMesh.isElementInMesh(idx, p, uv);
                    if (tf) {
                        std::vector<int> ids;
                        TopMesh.getElementIDs(idx, ids);
                        if (ids.size() == 3) {
                            top = TopNodes[ids[0]] * uv.x +
                                  TopNodes[ids[1]] * uv.y +
                                  TopNodes[ids[2]] * uv.z;
                            out = true;
                        } else if (ids.size() == 4) {
                            double N1, N2, N3, N4;
                            QuadShapeFunctions(uv.x, uv.y, N1, N2, N3, N4);
                            top = TopNodes[ids[0]] * N1 + TopNodes[ids[1]] * N2 +
                                  TopNodes[ids[2]] * N3 + TopNodes[ids[3]] * N4;
                            out = true;
                        }
                        break;
                    }
                }
            }
            if (bbot){
                CGAL::Circle_2<K> rc(cntr, BotRadius);
                BotSet.range_search(rc,std::back_inserter(LV));
                std::list<Vertex_handle2D>::const_iterator it = LV.begin();
                for (;it != LV.end(); ++it){
                    int idx = (*it)->info().id;
                    bool tf = TopMesh.isElementInMesh(idx, p, uv);
                    if (tf) {
                        std::vector<int> ids;
                        TopMesh.getElementIDs(idx, ids);
                        if (ids.size() == 3) {
                            bot = BotNodes[ids[0]] * uv.x +
                                  BotNodes[ids[1]] * uv.y +
                                  BotNodes[ids[2]] * uv.z;
                            out = true;
                        } else if (ids.size() == 4) {
                            double N1, N2, N3, N4;
                            QuadShapeFunctions(uv.x, uv.y, N1, N2, N3, N4);
                            bot = BotNodes[ids[0]] * N1 + BotNodes[ids[1]] * N2 +
                                  BotNodes[ids[2]] * N3 + BotNodes[ids[3]] * N4;
                            out = true;
                        }
                        break;
                    }
                }
            }
        }
        return out;
	}
}