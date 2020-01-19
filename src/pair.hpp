#ifndef PAIR_H_
#define PAIR_H_

#include <stdexcept>
#include <limits>
#include "vertex.h"
#include "vec4.hpp"
#include "mat4.hpp"
#include "map.hpp"
#include "face.hpp"

class Pair
{
private:
    Vec3 optPos;

public:
    double cost;
    int heapIndex;
    int index;
    int v[2];

    Pair()
    {
        cost = 0.0;
        optPos = Vec3();
        heapIndex = 0;
        index = 0;
        v[0] = v[1] = 0;
    }

    Pair(int v0, int v1)
    {
        cost = 0.0;
        optPos = Vec3();
        heapIndex = 0;
        index = 0;
        v[0] = v0;
        v[1] = v1;
    }

    Vec3 optimalPos() const
    {
        return optPos;
    }

    void updateOptiPos(const Vertex *vertices)
    {
        if (updateOptiPosInternal(vertices))
        {
            return;
        }

        // compute the (v0 + v1) / 2, v0 or v1
        Vec3 halfPos = (vertices[v[0]].p + vertices[v[1]].p) / 2;
        double cost_half = getCost(vertices, halfPos);
        double cost_v0 = getCost(vertices, vertices[v[0]].p);
        double cost_v1 = getCost(vertices, vertices[v[1]].p);

        optPos = halfPos;
        if (cost_v0 < cost_half)
        {
            optPos = vertices[v[0]].p;
        }
        if (cost_v1 < cost_half && cost_v1 < cost_v0)
        {
            optPos = vertices[v[1]].p;
        }
    }

    bool updateOptiPosInternal(const Vertex *vertices)
    {
        Mat4 A = vertices[v[0]].Q + vertices[v[1]].Q;
        A.set(3, 0, 0.0);
        A.set(3, 1, 0.0);
        A.set(3, 2, 0.0);
        A.set(3, 3, 1.0);
        Vec4 Y(0.0, 0.0, 0.0, 1.0);
        for (int i = 0; i < 4; ++i)
        {
            A.set(i, 3, -A.get(i, 3));
        }

        for (int i = 0; i < 3; ++i)
        {
            int j = 0;
            for (j = 0; j < 3; ++j)
            {
                if (std::abs(A.get(i, j)) >= 1e-6)
                {
                    break;
                }
            }
            if (j == 3)
                return false;
            for (int p = 0; p < 3; ++p)
            {
                if (p != i)
                {
                    double d = A.get(p, j) / A.get(i, j);
                    for (int k = 0; k < 4; ++k)
                    {
                        A.set(p, k, A.get(p, k) - A.get(i, k) * d);
                    }
                }
            }
        }

        for (int i = 0; i < 3; ++i)
        {
            int count = 0;
            for (int j = 0; j < 3; ++j)
            {
                if (std::abs(A.get(i, j)) < 1e-6)
                    count++;
            }
            if (count == 3)
                return false;
        }
        double index[3];
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (std::abs(A.get(i, j)) > 1e-6)
                {
                    index[j] = A.get(i, 3) / A.get(i, j);
                }
            }
        }
        optPos.x = index[0];
        optPos.y = index[1];
        optPos.z = index[2];
        return true;
    }

    void updateCost(const Vertex *vertices)
    {
        Vec4 y(optPos, 1.0);
        Mat4 A = vertices[v[0]].Q + vertices[v[1]].Q;
        Vec4 Ay = A * y;
        cost = y.dot(Ay);
    }

    double getFaceQuality(Vec3 v0, Vec3 v1, Vec3 v2)
    {
        Vec3 d10 = v1 - v0;
        Vec3 d20 = v2 - v0;
        Vec3 d12 = v1 - v2;
        Vec3 norm = d10.cross(d20);

        double len = norm.squareLength();
        if (len == 0.0)
        {
            return 0.0;
        }

        double b = d10.squareLength();
        if (b == 0.0)
        {
            return 0.0;
        }
        double t = b;
        t = d20.squareLength();
        if (b < t)
        {
            b = t;
        }

        t = d12.squareLength();
        if (b < t)
        {
            b = t;
        }
        return len / b;
    }

    void updateCost(Vertex *vertices, Map &faceMap)
    {
        // store old pair
        Vec3 v0_old = vertices[v[0]].p;
        Vec3 v1_old = vertices[v[1]].p;

        // normal
        std::vector<Vec3> onVec;

        // collect all normals
        // collect all faces in v0, skip faces with v1
        for (int i = 0; i < vertices[v[0]].neighbor.size(); ++i)
        {
            for (int j = i + 1; j < vertices[v[0]].neighbor.size(); ++j)
            {
                int neiIndex1 = vertices[v[0]].neighbor[i];
                int neiIndex2 = vertices[v[0]].neighbor[j];
                if (neiIndex1 == v[1] || neiIndex2 == v[1])
                {
                    continue;
                }
                Face realFace;
                int b = faceMap.get(Face(v[0], neiIndex1, neiIndex2), realFace);
                if (b)
                {
                    Vec3 norm = realFace.norm(vertices);
                    onVec.push_back(norm);
                }
            }
        }

        // collect all faces in v1, skip faces with v0
        for (int i = 0; i < vertices[v[1]].neighbor.size(); ++i)
        {
            for (int j = i + 1; j < vertices[v[1]].neighbor.size(); ++j)
            {
                int neiIndex1 = vertices[v[1]].neighbor[i];
                int neiIndex2 = vertices[v[1]].neighbor[j];
                if (neiIndex1 == v[0] || neiIndex2 == v[0])
                {
                    continue;
                }
                Face realFace;
                int b = faceMap.get(Face(v[1], neiIndex1, neiIndex2), realFace);
                if (b)
                {
                    Vec3 norm = realFace.norm(vertices);
                    onVec.push_back(norm);
                }
            }
        }

        // replace with the new position
        vertices[v[0]].p = optPos;
        vertices[v[1]].p = optPos;
        int i = 0;
        double MinCos = std::numeric_limits<double>::max();
        double MinQual = std::numeric_limits<double>::max();

        for (int i = 0; i < vertices[v[0]].neighbor.size(); ++i)
        {
            for (int j = i + 1; j < vertices[v[0]].neighbor.size(); ++j)
            {
                int neiIndex1 = vertices[v[0]].neighbor[i];
                int neiIndex2 = vertices[v[0]].neighbor[j];
                if (neiIndex1 == v[1] || neiIndex2 == v[1])
                {
                    continue;
                }
                Face currentFace;
                int b = faceMap.get(Face(v[0], neiIndex1, neiIndex2), currentFace);
                if (b)
                {
                    // get MinCos
                    Vec3 norm = currentFace.norm(vertices);
                    double ndiff = norm.dot(onVec[i++]);
                    MinCos = std::min(MinCos, ndiff);

                    // get face quality
                    double qt = getFaceQuality(vertices[v[0]].p, vertices[neiIndex1].p, vertices[neiIndex2].p);
                    MinQual = std::min(MinQual, qt);
                }
            }
        }

        for (int i = 0; i < vertices[v[1]].neighbor.size(); ++i)
        {
            for (int j = i + 1; j < vertices[v[1]].neighbor.size(); ++j)
            {
                int neiIndex1 = vertices[v[1]].neighbor[i];
                int neiIndex2 = vertices[v[1]].neighbor[j];
                if (neiIndex1 == v[1] || neiIndex2 == v[1])
                {
                    continue;
                }
                Face currentFace;
                int b = faceMap.get(Face(v[1], neiIndex1, neiIndex2), currentFace);
                if (b)
                {
                    // get MinCos
                    Vec3 norm = currentFace.norm(vertices);
                    double ndiff = norm.dot(onVec[i++]);
                    MinCos = std::min(MinCos, ndiff);

                    // get face quality
                    double qt = getFaceQuality(vertices[v[1]].p, vertices[neiIndex1].p, vertices[neiIndex2].p);
                    MinQual = std::min(MinQual, qt);
                }
            }
        }
    }

    double getCost(const Vertex *vertices, Vec3 pos)
    {
        Vec4 y(pos, 1.0);
        Mat4 A = vertices[v[0]].Q + vertices[v[1]].Q;
        Vec4 Ay = A * y;
        return y.dot(Ay);
    }

    friend bool operator<(const Pair &p1, const Pair &p2)
    {
        return p1.cost < p2.cost;
    }
};

#endif
