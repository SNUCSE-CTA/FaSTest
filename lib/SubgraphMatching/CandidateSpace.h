#pragma once
#include "SubgraphMatching/DataGraph.h"
#include "SubgraphMatching/PatternGraph.h"
#include "DataStructure/Graph.h"
#include "Base/Base.h"
#include "Base/BasicAlgorithms.h"

/**
 * @brief The Candidate Space structure
 * @date 2023-05
 */

namespace GraphLib {
namespace SubgraphMatching {
    BipartiteMaximumMatching BPSolver;
    enum STRUCTURE_FILTER {
        NO_STRUCTURE_FILTER,
        TRIANGLE_SAFETY,
        FOURCYCLE_SAFETY
    };
    enum NEIGHBOR_FILTER {
        NEIGHBOR_SAFETY,
        NEIGHBOR_BIPARTITE_SAFETY,
        EDGE_BIPARTITE_SAFETY
    };
    class SubgraphMatchingOption {
    public:
        STRUCTURE_FILTER structure_filter = FOURCYCLE_SAFETY;
        NEIGHBOR_FILTER neighborhood_filter = EDGE_BIPARTITE_SAFETY;
        int MAX_QUERY_VERTEX = 50, MAX_QUERY_EDGE = 250;
        long long max_num_matches = -1;
    };

    class CandidateSpace {
    public:
        SubgraphMatchingOption opt;
        CandidateSpace(DataGraph *data, SubgraphMatchingOption filter_option);

        ~CandidateSpace();

        CandidateSpace &operator=(const CandidateSpace &) = delete;

        CandidateSpace(const CandidateSpace &) = delete;

        inline int GetCandidateSetSize(const int u) const {return candidate_set_[u].size();};

        inline int GetCandidate(const int u, const int v_idx) const {return candidate_set_[u][v_idx];};

        bool BuildCS(PatternGraph *query);

        std::vector<int>& GetCandidates(int u) {
            return candidate_set_[u];
        }

        std::vector<int>& GetCandidateNeighbors(int cur, int cand_idx, int nxt) {
            return candidate_neighbors[cur][cand_idx][query_->GetAdjIdx(cur, nxt)];
        }
        int GetCandidateNeighbor(int cur, int cand_idx, int nxt, int nxt_idx) {
            return candidate_neighbors[cur][cand_idx][query_->GetAdjIdx(cur, nxt)][nxt_idx];
        }

        dict GetCSInfo() {return CSInfo;};

        int GetNumCSVertex() {return num_candidate_vertex;};
        int GetNumCSEdge() {return num_candidate_edge;};
    private:
        dict CSInfo;
        DataGraph *data_;
        PatternGraph *query_;
        std::vector<std::vector<std::vector<std::vector<int>>>> candidate_neighbors;

        std::vector<std::vector<int>> candidate_set_;
        std::vector<int> neighbor_label_frequency;
        int num_candidate_vertex = 0, num_candidate_edge = 0;
        bool* in_neighbor_cs;
        bool** BitsetCS;
        bool** BitsetEdgeCS;
        int* num_visit_cs_;

        bool BuildInitialCS();

        void ConstructCS();

        bool InitRootCandidates(int root);

        bool RefineCS();

        void PrepareNeighborSafety(int cur);

        bool CheckNeighborSafety(int cur, int cand);

        bool NeighborBipartiteSafety(int cur, int cand);

        bool EdgeBipartiteSafety(int cur, int cand);

        inline bool EdgeCandidacy(int query_edge_id, int data_edge_id);

        bool TriangleSafety(int query_edge_id, int data_edge_id);

        bool FourCycleSafety(int query_edge_id, int data_edge_id);

        bool CheckSubStructures(int cur, int cand);

        bool NeighborFilter(int cur, int cand) {
            switch (opt.neighborhood_filter) {
                case NEIGHBOR_SAFETY:
                    return CheckNeighborSafety(cur, cand);
                case NEIGHBOR_BIPARTITE_SAFETY:
                    return NeighborBipartiteSafety(cur, cand);
                case EDGE_BIPARTITE_SAFETY:
                    return EdgeBipartiteSafety(cur, cand);
            }
        };

        bool StructureFilter(int query_edge_id, int data_edge_id) {
            switch (opt.structure_filter) {
                case NO_STRUCTURE_FILTER:
                    return true;
                case TRIANGLE_SAFETY:
                    return TriangleSafety(query_edge_id, data_edge_id);
                case FOURCYCLE_SAFETY:
                    return TriangleSafety(query_edge_id, data_edge_id) and FourCycleSafety(query_edge_id, data_edge_id);
            }
            return true;
        };

    };

    CandidateSpace::CandidateSpace(DataGraph *data, SubgraphMatchingOption filter_option) {
        opt = filter_option;
        data_ = data;
        BitsetCS = new bool*[opt.MAX_QUERY_VERTEX];
        for (int i = 0; i < opt.MAX_QUERY_VERTEX; i++) {
            BitsetCS[i] = new bool[data->GetNumVertices()]();
        }
        BitsetEdgeCS = new bool*[opt.MAX_QUERY_EDGE];
        for (int i = 0; i < opt.MAX_QUERY_EDGE; i++) {
            BitsetEdgeCS[i] = new bool[data->GetNumEdges()]();
        }
        in_neighbor_cs = new bool[data->GetNumVertices()]();
        neighbor_label_frequency.resize(data->GetNumVertices());
        num_visit_cs_ = new int[data_->GetNumVertices()]();
        candidate_set_.resize(opt.MAX_QUERY_VERTEX);
        BPSolver.Initialize(50, data_->GetMaxDegree(), 50);
        fprintf(stderr, "Constructing Candidate Space: %d %d\n", opt.MAX_QUERY_VERTEX, opt.MAX_QUERY_EDGE);
    }

    CandidateSpace::~CandidateSpace() {
        for (int i = 0; i < opt.MAX_QUERY_VERTEX; i++) {
            delete[] BitsetCS[i];
        }
        delete[] BitsetCS;
        for (int i = 0; i < opt.MAX_QUERY_EDGE; i++) {
            delete[] BitsetEdgeCS[i];
        }
        delete[] BitsetEdgeCS;
        delete[] num_visit_cs_;
        delete[] in_neighbor_cs;
    }


    bool CandidateSpace::BuildCS(PatternGraph *query) {
        CSInfo.clear();
        query_ = query;
        num_candidate_vertex = num_candidate_edge = 0;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            memset(BitsetCS[i], false, data_->GetNumVertices());
        }
        for (int i = 0; i < query_->GetNumEdges(); i++) {
            memset(BitsetEdgeCS[i], false, data_->GetNumEdges());
        }
        memset(num_visit_cs_, 0, data_->GetNumVertices());
        BPSolver.Reset();
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            candidate_set_[i].clear();
        }
        Timer CStimer; CStimer.Start();
        BuildInitialCS();
        RefineCS();
        ConstructCS();
        CStimer.Stop();
        CSInfo["CSBuildTime"] = CStimer.GetTime();
        CSInfo["#CSVertex"] = num_candidate_vertex;
        CSInfo["#CSEdge"] = num_candidate_edge;
        return true;
    }

    void CandidateSpace::ConstructCS() {
        candidate_neighbors.clear();
        candidate_neighbors.resize(query_->GetNumVertices());
        std::vector <int> candidate_index(data_->GetNumVertices());
        for (int i = 0; i < query_->GetNumVertices(); ++i) {
            candidate_neighbors[i].resize(GetCandidateSetSize(i));
        }
        for (int u = 0; u < query_->GetNumVertices(); u++) {
            int u_label = query_->GetVertexLabel(u);
            int u_degree = query_->GetDegree(u);
            num_candidate_vertex += GetCandidateSetSize(u);
            for (int idx = 0; idx < GetCandidateSetSize(u); idx++) {
                candidate_index[candidate_set_[u][idx]] = idx;
                candidate_neighbors[u][idx].resize(query_->GetDegree(u));
            }

            for (int uc : query_->GetNeighbors(u)) {
                int query_edge_idx = query_->GetEdgeIndex(uc, u);
                for (int vc_idx = 0; vc_idx < candidate_set_[uc].size(); ++vc_idx) {
                    int vc = candidate_set_[uc][vc_idx];
                    for (int data_edge_idx : data_->GetIncidentEdges(vc, u_label)) {
                        int v = data_->GetOppositePoint(data_edge_idx);
                        if (data_->GetDegree(v) < u_degree) break;
                        if (!BitsetEdgeCS[query_edge_idx][data_edge_idx]) continue;
                        num_candidate_edge++;
                        candidate_neighbors[u][candidate_index[v]][query_->GetAdjIdx(u, uc)].emplace_back(vc_idx);
                    }
                }
            }
        }
    }

} }
