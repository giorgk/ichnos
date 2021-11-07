//
// Created by giorg on 11/4/2021.
//

#include "ichnos_utils.h"

#ifndef ICHNOS_ICHNOS_MESH_H
#define ICHNOS_ICHNOS_MESH_H

namespace ICHNOS{
    class Mesh2D {
    public:
        Mesh2D(){};
        void setNodes(std::vector<vec3> &nd);
        void setMesh(std::vector<std::vector<int>> &msh, bool zerobased);
        void CalculateCentroid(std::vector<vec3> &cc);
        void CalculateCentroid(int elid, vec3 &cc);
        bool isElementInMesh(int elid, vec3 &p, vec3 &uv);
        void getElementIDs(int elid, std::vector<int> &ids);
        vec3& getNode(int el, int idx);
        void setInverseTransformTolerance(double t){invtol = t;}

    private:
        std::vector<vec3> nodes;
        std::vector<std::vector<int>> mesh;
        double invtol = 0.000001;
    };

    void Mesh2D::setNodes(std::vector<vec3> &nd) {
        nodes = nd;
    }

    void Mesh2D::setMesh(std::vector<std::vector<int>> &msh, bool zerobased) {
        for (unsigned int i = 0; i < msh.size(); ++i){
            std::vector<int> tmp;
            for (unsigned int j = 0; j < msh[i].size(); ++j){
                if (zerobased){
                    if (msh[i][j] >= 0)
                        tmp.push_back(msh[i][j]);
                }
                else{
                    if (msh[i][j] > 0)
                        tmp.push_back(msh[i][j] - 1);
                }
            }
            mesh.push_back(tmp);
        }
    }

    void Mesh2D::CalculateCentroid(std::vector<vec3> &cc) {
        cc.clear();
        for (unsigned int i = 0; i < mesh.size(); ++i){
            vec3 bc;
            CalculateCentroid(i, bc);
            cc.push_back(bc);
        }
    }

    void Mesh2D::CalculateCentroid(int elid, vec3 &cc) {
        cc.zero();
        for (unsigned int i = 0; i < mesh[elid].size(); ++i){
            cc.x += nodes[mesh[elid][i]].x;
            cc.y += nodes[mesh[elid][i]].y;
            cc.z += nodes[mesh[elid][i]].z;
        }
        double n = static_cast<double>(mesh[elid].size());
        cc.x = cc.x/n;
        cc.y = cc.y/n;
        cc.z = cc.z/n;
    }

    bool Mesh2D::isElementInMesh(int elid, vec3 &p, vec3 &uv) {
        if (elid >= mesh.size())
            return false;

        if (mesh[elid].size() == 3){
            return isPointInTriangle(p,
                                     nodes[mesh[elid][0]],
                                     nodes[mesh[elid][1]],
                                     nodes[mesh[elid][2]],
                                     uv);

        }
        else if (mesh[elid].size() == 4){
            bool tf1, tf2;
            tf1 = isPointInTriangle(p,
                                     nodes[mesh[elid][0]],
                                     nodes[mesh[elid][1]],
                                     nodes[mesh[elid][2]],
                                     uv);
            if (!tf1){
                tf2 = isPointInTriangle(p,
                                        nodes[mesh[elid][0]],
                                        nodes[mesh[elid][2]],
                                        nodes[mesh[elid][3]],
                                        uv);
            }
            if (!tf1 && !tf2){
                return false;
            }
            else{
                return isPointInQuad(p,
                                     nodes[mesh[elid][0]],
                                     nodes[mesh[elid][1]],
                                     nodes[mesh[elid][2]],
                                     nodes[mesh[elid][3]],
                                     uv,invtol);
            }
        }
        else{
            return false;
        }
        return false;
    }

    void Mesh2D::getElementIDs(int elid, std::vector<int> &ids) {
        ids.clear();
        if (elid >= mesh.size())
            return;
        else
            ids = mesh[elid];
    }
}

#endif //ICHNOS_ICHNOS_MESH_H
