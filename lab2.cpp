#include <iostream>
#include <list>
#include <map>
#include <set>
#include <unordered_set>

#include <glm/gtc/matrix_inverse.hpp>
#include <spdlog/spdlog.h>

#include "Labs/2-GeometryProcessing/DCEL.hpp"
#include "Labs/2-GeometryProcessing/tasks.h"

namespace VCX::Labs::GeometryProcessing {

#include "Labs/2-GeometryProcessing/marching_cubes_table.h"

    /******************* 1. Mesh Subdivision *****************/
    void SubdivisionMesh(Engine::SurfaceMesh const & input, Engine::SurfaceMesh & output, std::uint32_t numIterations) {
        // your code here
        output = Engine::SurfaceMesh(input);
        Engine::SurfaceMesh temp = Engine::SurfaceMesh(input);
        for (int i = 0; i < numIterations; ++i){
            output.Positions.clear();
            output.Indices.clear();
            int  vercount = temp.Positions.size(), facecount = temp.Indices.size() / 3;
            DCEL links = DCEL(vercount, facecount);
            (&links)->AddFaces(temp.Indices);
            int edgecount = links.GetEdges().size();
            for (int i = 0; i < vercount; ++i){
                DCEL::Vertex v = links.GetVertex(i);
                auto neighbors = v.GetNeighbors();
                int n = neighbors.size();
                float u = 3.0 / 16;
                if (n != 3)
                    u = 3.0 / (8 * n);
                glm::vec3 pos = (1 - n * u) * temp.Positions[i];
                for (auto id : neighbors) {
                    glm::vec3 c = temp.Positions[id];
                    pos += u * c;
                }
                output.Positions.push_back(pos);
            }
            std::map<std::pair<int, int>, int> tmp;
            for (auto p : links.GetEdges()) {
                int u = p->From(), v = p->To();
                tmp[std::make_pair(u, v)] = tmp[std::make_pair(v, u)] = vercount++;
                float c1 = 3.0 / 8, c2 = 1.0 / 8;
                glm::vec3 pos = c1 * (temp.Positions[u] + temp.Positions[v]) + c2 * (temp.Positions[p->OppositeVertex()] + temp.Positions[p->PairEdge()->OppositeVertex()]);
                output.Positions.push_back(pos);
            }
            int A[3];
            for (DCEL::Triangle const & f : links.GetFaces()) {
                for (int i = 0; i <= 2; ++i) {
                    A[i] = tmp[std::make_pair(f.Edges(i)->From(), f.Edges(i)->To())];
                }
                for (int i = 0; i <= 2; ++i) {
                    output.Indices.push_back(A[i]);
                    output.Indices.push_back(f.Edges(i)->To());
                    output.Indices.push_back(A[(i + 1) % 3]);
                }
                for (int i = 0; i <= 2; ++i) {
                    output.Indices.push_back(A[i]);
                }
            }
            temp = output;
        }
    }

    /******************* 2. Mesh Parameterization *****************/
    void Parameterization(Engine::SurfaceMesh const & input, Engine::SurfaceMesh & output, const std::uint32_t numIterations) {
        // your code here
        DCEL links;
        links.AddFaces(input.Indices);
        int bgn = 0, n = input.Positions.size();
        for (int i = 0; i < n; ++i){
            DCEL::Vertex v = links.GetVertex(i);
            if (v.IsSide()){
                bgn = i;
                break;
            }
        }
        int lst = -1, cur = bgn;
        std::vector<int> border;
        border.push_back(bgn);
        while (lst == -1 || cur != bgn){
            DCEL::Vertex v = links.GetVertex(cur);
            if (lst == -1 || v.GetSideNeighbors().first == lst){
                lst = cur;
                cur = v.GetSideNeighbors().second;
            }
            else{
                lst = cur;
                cur = v.GetSideNeighbors().first;
            }
            border.push_back(cur);
        }
        int bsize = border.size();
        float ac = 2 * acos(-1) / bsize;
        std::vector<glm::vec2> uv;
        for (int i = 0; i < n; ++i)
            uv.push_back(glm::vec2(0,0));
        for (int i = 0; i < bsize; ++i){
            uv[border[i]].x = 0.5 + 0.5 * cos(i * ac);
            uv[border[i]].y = 0.5 + 0.5 * sin(i * ac);
        }
        for (int k = 0; k < numIterations; ++k){
            std::vector<glm::vec2> tmp;
            for (int i = 0; i < n; ++i)
                tmp.push_back(uv[i]);
            for (int i = 0; i < n; ++i){
                DCEL::Vertex v = links.GetVertex(i);
                if (v.IsSide())
                    continue;
                int d = v.GetNeighbors().size();
                float lam = 1.0 / d;
                uv[i] = glm::vec2(0,0);
                for (int j = 0; j < d; ++j)
                    uv[i] += tmp[v.GetNeighbors()[j]];
                uv[i] *= lam;
            }
        }
        for (int i = 0; i < n; ++i){
            output.TexCoords.push_back(uv[i]);
            output.Positions.push_back(input.Positions[i]);
        }
        for (int i = 0; i < input.Indices.size(); ++i)
            output.Indices.push_back(input.Indices[i]);
    }


    /******************* 3. Mesh Simplification *****************/
    struct node {
        int i, j;
        float c;
        glm::vec3 v;
        node(int _i, int _j, float _c, glm::vec3 _v) { i = _i, j = _j, c = _c, v = _v; }
    };
    glm::mat4 Q[52000];
    node func(Engine::SurfaceMesh const & input, int i, int j) {
        glm::mat4 Qp = Q[i] + Q[j];
        glm::mat4 Qu = {
            Qp[0], Qp[1], Qp[2], {0, 0, 0, 1}
        };
        Qu = glm::transpose(Qu);
        glm::vec4 vu;
        if (glm::determinant(Qu) < 1e-6 || 1) {
            glm::vec3 vs = (input.Positions[i] + input.Positions[j]) / (float) 2.0;
            vu = { vs[0], vs[1], vs[2], 1 };
        } 
        else {
            glm::mat4 Qi = glm::inverse(Qu);
            vu = { Qi[3][0], Qi[3][1], Qi[3][2], 1 };
        }
        glm::vec3 vs = (input.Positions[i] + input.Positions[j]) / (float) 2.0;
        float vl = glm::dot(vu, Qp * vu);
        glm::vec3 vk = { vu[0], vu[1], vu[2] };
        return node(i, j, vl, vk);
    }
    int ft[52000];
    int getf(int k) {
        if (ft[k] == 0 || ft[k] == k) return k;
        return ft[k] = getf(ft[k]);
    }
    std::map<std::pair<int, std::pair<int, int>>, bool> bz;
    std::list<node> wlist;
    int index[52000];

    void SimplifyMesh(Engine::SurfaceMesh const & input, Engine::SurfaceMesh & output, float valid_pair_threshold, float simplification_ratio) {
        // your code here
        memset(ft, 0, sizeof(ft));
        int  vercount = input.Positions.size(), facecount = input.Indices.size() / 3;
        DCEL links = DCEL(vercount, facecount);
        (&links)->AddFaces(input.Indices);
        for (int i = 0; i < vercount; ++i) {
            Q[i] = glm::mat4(0);
        }
        for (int i = 0; i < facecount; ++i) {
            int v[3] = { input.Indices[3 * i], input.Indices[3 * i + 1], input.Indices[3 * i + 2] };
            glm::vec3 v1 = input.Positions[v[1]] - input.Positions[v[0]];
            glm::vec3 v2 = input.Positions[v[2]] - input.Positions[v[1]];
            glm::vec3 norm = glm::cross(v1, v2);
            norm = norm / float(sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]));
            float d = -glm::dot(norm, input.Positions[v[0]]);
            glm::mat4 p = {{norm[0], norm[1], norm[2], d}, { 0, 0, 0, 0}, { 0, 0, 0, 0}, { 0, 0, 0, 0}};
            glm::mat4 kp = p * glm::transpose(p);
            for (int j = 0; j <=3; ++j)
                Q[v[j]] += kp;
        }
        wlist.clear();
        for (int i = 0; i < vercount; ++i)
            for (int j = i+1; j < vercount; ++j) {
            bool bz = 0;
            DCEL::Vertex ei = links.GetVertex(i);
            for (auto ej : ei.GetNeighbors()) {
                if (ej == j) {
                    bz = 1;
                    break;
                }
            }
            glm::vec3 tmp = input.Positions[i] - input.Positions[j];
            if (! bz && sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]) > valid_pair_threshold)
                continue;
            node nd = func(input, i, j);
            wlist.push_back(nd);
        }
        int l = input.Positions.size() * (1 - simplification_ratio);
        Engine::SurfaceMesh temp = input;
        for (int i = 0; i < l; ++i) {
            float minn = 1800016527;
            int mi, mj;
            glm::vec3 mu;
            int cnt = 0;
            for (std::list<node>::iterator iter = wlist.begin(); iter != wlist.end();) {
                int p = iter->i, q = iter->j;
                if (p != getf(p) || q != getf(q)) {
                    p = getf(p), q = getf(q);
                    if (p == q) {
                        iter = wlist.erase(iter);
                        continue;
                    }
                    *iter = func(temp, p, q);
                }
                if (iter->c < minn) {
                    minn = iter->c;
                    mi = p, mj = q;
                    mu = iter->v;
                }
                iter++;
            }
            if (mi > mj) std::swap(mi, mj);
            ft[mj] = mi;
            temp.Positions[mi] = mu;
            Q[mi] += Q[mj];
        }
        output.Positions.clear();
        output.Indices.clear();
        for (int i = 0; i < vercount; ++i) {
            if (getf(i) == i) {
                output.Positions.push_back(temp.Positions[i]);
            }
        }
        int cnt = 0;
        for (int i = 0; i < facecount; ++i) {
            index[i] = cnt;
            if (getf(i) == i) cnt++;
        }
        bz.clear();
        for (int i = 0; i < facecount; ++i) {
            int v[3] = { getf(input.Indices[3 * i]), getf(input.Indices[3 * i + 1]), getf(input.Indices[3 * i + 2]) };
            if (v[0] != v[1] && v[1] != v[2] && v[0] != v[2]) {
                for (int j = 0; j < 3; ++j)
                    output.Indices.push_back(index[v[j]]);
            }
        }
    }

    /******************* 4. Mesh Smoothing *****************/
    double sum_square(glm::vec3 x, glm::vec3 y) {
        return x.x * y.x + x.y * y.y + x.z * y.z;
    }
    double sum_len(glm::vec3 x){
        return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    }
    float func(glm::vec3 x, glm::vec3 y, glm::vec3 p, glm::vec3 q) {
        float alpha = sum_square((x - p), (y - p)) / (sum_len(x - p) * sum_len(y - p));
        float beta = sum_square((x - q), (y - q)) / (sum_len(x - q) * sum_len(y - q));
        return std::max(std::min(alpha / std::sqrt(1 - alpha * alpha) + beta / std::sqrt(1 - beta * beta),800.0f),0.0f);
    }
    void SmoothMesh(Engine::SurfaceMesh const & input, Engine::SurfaceMesh & output, std::uint32_t numIterations, float lambda, bool useUniformWeight) {
        // your code here
        std::vector<glm::vec3> P;
        std::vector<glm::vec3> P1;
        std::vector<float> wsum;
        int n = input.GetVertexCount();
        for (int i = 0; i < n; ++i) {
            P.push_back(input.Positions[i]);
            P1.push_back(glm::vec3(0,0,0));
            wsum.push_back(0);
        }
        DCEL links;
        links.AddFaces(input.Indices);
        for (int k = 0; k < numIterations; ++k){
            for (int i = 0; i < n; ++i){
                wsum[i] = 0;
                P1[i] = glm::vec3(0, 0, 0);
            }
            for (DCEL::HalfEdge const * e : links.GetEdges()) {
                float wij = 1;
                DCEL::Vertex pu = links.GetVertex(e->From()), pv = links.GetVertex(e->To());
                if (useUniformWeight==0) 
                    wij = func(P[e->From()], P[e->To()], P[e->OppositeVertex()], P[e->PairOppositeVertex()]);
                wsum[e->From()] += wij;
                wsum[e->To()] += wij;
                glm::vec3 u = P[e->To()], v = P[e->From()];
                u *= wij;
                v *= wij;
                P1[e->From()] += u;
                P1[e->To()] += v;
            }
            for (int i = 0; i < n; ++i){
                P[i] *= (1 - lambda);
                P1[i] /= wsum[i];
                P1[i] *= lambda;
                P[i] += P1[i];
            }
        }
        for (int i = 0; i < n; ++i)
            output.Positions.push_back(P[i]);
        for (int i = 0; i < input.Indices.size(); ++i)
            output.Indices.push_back(input.Indices[i]);
    }

    /******************* 5. Marching Cubes *****************/
    glm::vec3 get_pos(glm::vec3 u,float valu,float valv,glm::vec3 unit) {
        float sum = valu - valv;
        glm::vec3 res = u + (valu / sum) * unit;
        return res;
    }
    bool check(glm::vec3 p, double f) {
        return std::fabs(p.x) < f / 1000 && std::fabs(p.y) < f / 1000 && std::fabs(p.z) < f / 1000;
    }
    void MarchingCubes(Engine::SurfaceMesh & output, const std::function<float(const glm::vec3 &)> & sdf, const glm::vec3 & grid_min, const float dx, const int n) {
        // your code here
        std::vector<int> * Pointid;
        Pointid = new std::vector<int>[n * n * n];
        glm::vec3 unit[3];
        unit[0] = glm::vec3(1, 0, 0);
        unit[1] = glm::vec3(0, 1, 0);
        unit[2] = glm::vec3(0, 0, 1);
        int id[12];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < n; y++)
                for (int z = 0; z < n; z++){
                    int cubeid = x * n * n + y * n + z;
                    Pointid[cubeid].clear();
                    glm::vec3 v[8];
                    for (int i = 0; i < 8; i++)
                        v[i] = glm::vec3(grid_min.x + 1.0*x * dx + (i & 1) * dx, grid_min.y + 1.0*y * dx + (i >> 1 & 1) * dx, grid_min.z + 1.0*z * dx + (i >> 2 & 1) * dx);
                    uint32_t state = 0;
                    for (int i = 0; i < 8; i++)
                        if (sdf(v[i]) <= 0)
                            state |= (1 << i);
                    uint32_t Estate = c_EdgeStateTable[state];
                    int tot = 0;
                    for (int i = 0; i < 12; ++i){
                        glm::vec3 sta = v[0] + dx * (i & 1) * unit[((i >> 2) + 1 ) % 3] + dx * (i >> 1 & 1) * unit[((i >> 2) + 2) % 3];
                        glm::vec3 end = sta + dx * unit[i >> 2];
                        glm::vec3 newp = get_pos(sta, sdf(sta), sdf(end), dx * unit[i >> 2]);
                        int cid;
                        int pid = -1;
                        if (x != 0){
                            cid = (x - 1) * n * n + y * n + z;
                                for (int j = 0; j < Pointid[cid].size(); j++)
                                    if (check(newp - output.Positions[Pointid[cid][j]], dx)) {
                                        pid = Pointid[cid][j];
                                        break;
                                    }
                        }
                        if (y != 0) {
                                cid = x * n * n + (y - 1) * n + z;
                                for (int j = 0; j < Pointid[cid].size(); j++)
                                    if (check(newp - output.Positions[Pointid[cid][j]], dx)) {
                                        pid = Pointid[cid][j];
                                        break;
                                    }
                            }
                        if (z != 0) {
                                cid = x * n * n + y * n + z - 1;
                                for (int j = 0; j < Pointid[cid].size(); j++)
                                    if (check(newp - output.Positions[Pointid[cid][j]],dx)) {
                                        pid = Pointid[cid][j];
                                        break;
                                    }
                            }
                            if (pid == -1){
                                pid=output.Positions.size();
                                output.Positions.push_back(newp);
                            }
                            Pointid[cubeid].push_back(pid);
                            id[i] = pid;
                    }
                    for (int j = 0; j < 12; j++) {
                        if (c_EdgeOrdsTable[state][j] != -1) {
                            output.Indices.push_back(id[c_EdgeOrdsTable[state][j]]);
                        }
                    }
                }
    }
} // namespace VCX::Labs::GeometryProcessing