#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <set>
#include <memory.h>

#include "mesh.h"

Mesh::Mesh() {
    vertices = new Vertex[MAX_VERTEX];
    pairs = new Pair[MAX_PAIR];
    vOffset = 0;
    pOffset = 0;
    faceCount = 0;
    vertexCount = 0;    
    t = 0;
    rate = 1.0;
    valid = new bool[MAX_VERTEX];
    inPair = new bool[MAX_VERTEX];
    inFace = new bool[MAX_VERTEX];
    memset(valid, 0, MAX_VERTEX * sizeof(bool));
    memset(inPair, 0, MAX_VERTEX * sizeof(bool));
    memset(inFace, 0, MAX_VERTEX * sizeof(bool));
}

Mesh::~Mesh() {
    if (valid) delete[] valid;
    if (inPair) delete[] inPair;
    if (inFace) delete[] inFace;
    if (vertices) delete[] vertices;
    if (pairs) delete[] pairs;
}

void Mesh::load(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) {
        std::cout << "Not such file" << std::endl;
        exit(1);
    }
    char line[1000];
    char type;
    while (fin >> line) {
        type = line[0];
        if (type == 'v') {
            double x = 0, y = 0, z = 0;
            sscanf(line, "%c%lf%lf%lf", &type, &x, &y, &z);
            Vec3 p(x, y, z);
            addVertex(p);
        }
        else if (type == 'f'){
            int indices[3];
            sscanf(line, "%c%d%d%d", &type, &indices[0], &indices[1], &indices[2]);
            Face f(indices);
            addFace(f);
        }
        else {
            //type == '#'
        }
    }
}

void Mesh::save(const std::string& path) {
    std::ofstream fout(path);
    int vNum = 0, fNum = 0;
    for (int index = 0; index < vOffset; ++index) {
        if (valid[index]) {
            fout << "v " << vertices[index].p.x << " " << vertices[index].p.y << " " << vertices[index].p.z << std::endl;
            ++vNum;
        }
    }
    assert(vNum == vertexCount);
    for (int index = 0; index < vOffset; ++index) {
        if (valid[index]) {
            for (int i = 0; i < vertices[index].neighbor.size(); ++i) {
                int neiIndex1 = vertices[index].neighbor[i];
                for (int j = i + 1; j < vertices[index].neighbor.size(); ++j) {
                    int neiIndex2 = vertices[index].neighbor[j];
                    if (!inFace[neiIndex1] && !inFace[neiIndex2]) {
                        if (vertices[neiIndex1].isNeighbor(neiIndex2)) {
                            ++fNum;
                            fout << "f " << index << " " << neiIndex1 << " " << neiIndex2 << std::endl;
                        }
                    }
                }
            }
            inFace[index] = true;
        }
    }
    fout.seekp(0, std::ios::beg);
    fout << "# " << vNum << " vertices, " << fNum << " faces" << std::endl;
}

int Mesh::addVertex(const Vec3& p) {
    int index = vOffset;
    vertices[vOffset].setPos(p);
    valid[vOffset] = true;
    ++vOffset;
    ++vertexCount;
    return index;
}

void Mesh::delVertex(int index) {
    valid[index] = false;
    --vertexCount;
}

void Mesh::addFace(const Face& f) {
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (!vertices[f.indices[i]].isNeighbor(f.indices[j])) {
                vertices[f.indices[i]].addNeighbor(f.indices[j]);
                vertices[f.indices[j]].addNeighbor(f.indices[i]);
            }
        }
    }
    ++faceCount;
}

int Mesh::addPair(int v1, int v2) {
    int pairIndex = pOffset;
    pairs[pOffset].v[0] = v1;
    pairs[pOffset].v[1] = v2;
    pairs[pOffset].index = pairIndex;
    ++pOffset;
    return pairIndex;
}

void Mesh::computeQ() {
    for (int index = 0; index < vertexCount; ++index) {
        vertices[index].computeQ(vertices);
    }
}

void Mesh::computeValidPairs() {
    for (int index = 0; index < vertexCount; ++index) {
        for (int i = 0; i < vertices[index].neighbor.size(); ++i) {
            int neighborIndex = vertices[index].neighbor[i];
            if (!inPair[neighborIndex]) {
                int pairIndex = addPair(index, neighborIndex);
                pairs[pairIndex].updateOptiPos(vertices);
                pairs[pairIndex].updateCost(vertices);
            }
        }
        // there may also be |v1 - v2| < t
        inPair[index] = true;
    }
}

void Mesh::simplify() {
    computeQ();
    std::cout << "[info] Compute Q end" << std::endl;
    computeValidPairs();
    std::cout << "[info] Compute valid pairs end" << std::endl;
    heap.build(pairs, pOffset);
    std::cout << "[info] Heap build end" << std::endl;

    int nowCount = faceCount;
    int iter = 0;
    while (nowCount > int(rate * faceCount)) {
        if (iter % 100 == 0) {
            std::cout << "[info] Iter: " << iter << " count: " << nowCount << std::endl;
        }
        int pairIndex = heap.top();
        if (pairIndex < 0) {
            break;
        }
        Pair minPair = pairs[heap.top()];
        heap.del();
        update(minPair);
        nowCount -= 2;
        iter++;
    }
}

void Mesh::update(const Pair& pair) {
    //optimize: use v[0] to store new vertex
    Vec3 newPos = pair.optimalPos();
    int newIndex = pair.v[0];
    vertices[newIndex].setPos(newPos);
    vertices[newIndex].Q += vertices[pair.v[1]].Q;

    /* get new neighbor */
    //std::set<int> newNeighbor;
    //for (int i = 0; i < vertices[pair.v[0]].neighbor.size(); ++i) {
    //    int neighborIndex = vertices[pair.v[0]].neighbor[i];
    //    if (neighborIndex != pair.v[1]) {
    //        newNeighbor.insert(neighborIndex);
    //    }
    //}
    for (int i = 0; i < vertices[pair.v[1]].neighbor.size(); ++i) {
        int neighborIndex = vertices[pair.v[1]].neighbor[i];
        if (neighborIndex != pair.v[0]) {
            if (!vertices[newIndex].isNeighbor(neighborIndex)) {
                vertices[newIndex].addNeighbor(neighborIndex);
                vertices[neighborIndex].addNeighbor(newIndex);
            }
            vertices[neighborIndex].delNeighbor(pair.v[1]);
        }
    }
    vertices[newIndex].delNeighbor(pair.v[1]);
    //for (int i = 0; i < vertices[newIndex].neighbor.size(); ++i) {
    //    int neighborIndex = vertices[newIndex].neighbor[i];
    //    vertices[neighborIndex].delNeighbor(pair.v[1]);
    //    if (vertices[neighborIndex].isNeighbor)
    //    vertices[neighborIndex].addNeighbor(newIndex);
    //}
    //std::set<int>::iterator iter;
    ////for (iter = newNeighbor.begin(); iter != newNeighbor.end(); ++iter) {
    ////    vertices[newIndex].addNeighbor(*iter);
    ////}
    //for (iter = newNeighbor.begin(); iter != newNeighbor.end(); ++iter) {
    //    vertices[*iter].delNeighbor(pair.v[0]);
    //    vertices[*iter].delNeighbor(pair.v[1]);
    //    vertices[*iter].addNeighbor(newIndex);
    //    //vertices[*iter].computeQ(vertices);
    //}

    /*update pairs*/    
    //add v[1] pairs to new vertex(v[0])
    for (int i = 0; i < vertices[pair.v[1]].pairs.size(); ++i) {
        int pairIndex = vertices[pair.v[1]].pairs[i];
        if (pairs[pairIndex].v[0] == pair.v[1]) {
            if (pairs[pairIndex].v[1] == pair.v[0]) {
                //pair between v[0] and v[1]
                assert(pairIndex == pair.index);
                vertices[newIndex].delPair(pairIndex);
                heap.remove(pairs[pairIndex]);
                continue;
            }
            else {
                pairs[pairIndex].v[0] = newIndex;
            }
        }
        else {
            if (pairs[pairIndex].v[0] == pair.v[0]) {
                //pair between v[0] and v[1]
                assert(pairIndex == pair.index);
                vertices[newIndex].delPair(pairIndex);
                heap.remove(pairs[pairIndex]);
                continue;
            }
            else {
                pairs[pairIndex].v[1] = newIndex;
            }
        }
        if (vertices[newIndex].hasPair(pairIndex, pairs)) {
            heap.remove(pairs[pairIndex]);
        }
        else {
            vertices[newIndex].addPair(pairIndex);
        }
    }
    //update cost & optimal pos
    for (int i = 0; i < vertices[newIndex].pairs.size(); ++i) {
        int pairIndex = vertices[newIndex].pairs[i];
        pairs[pairIndex].updateOptiPos(vertices);
        pairs[pairIndex].updateCost(vertices);
        heap.update(pairs[pairIndex]);
    }

    //heap.remove(pair);
    delVertex(pair.v[1]);
}