#include "SubgraphMatching/CandidateSpace.h"

namespace GraphLib::SubgraphMatching {
    bool CandidateSpace::BuildInitialCS() {
        std::memset(num_visit_cs_, 0, data_->GetNumVertices() * sizeof(int));
        std::vector<int> initial_cs_size(query_->GetNumVertices(), 0);
        std::vector<std::vector<int>> built_neighbors(query_->GetNumVertices());
        int root = 0;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            int l = query_->GetVertexLabel(i);
            int d = query_->GetDegree(i);
            initial_cs_size[i] = data_->GetVerticesByLabel(l).size();
            if (initial_cs_size[i] <= initial_cs_size[0]) {
                root = i;
            }
        }
        InitRootCandidates(root);
        for (int uc : query_->GetNeighbors(root)) {
            built_neighbors[uc].push_back(root);
        }
        for (int i = 1; i < query_->GetNumVertices(); i++) {
            int cur = -1;
            for (int j = 0; j < query_->GetNumVertices(); j++) {
                if (!candidate_set_[j].empty()) continue;
                if (cur == -1) { cur = j; continue; }
                if (built_neighbors[j].size() > built_neighbors[cur].size()) { cur = j; continue; }
                if (built_neighbors[j].size() == built_neighbors[cur].size() and
                    initial_cs_size[j] < initial_cs_size[cur]) { cur = j; continue; }
            }
            int cur_label = query_->GetVertexLabel(cur);
            int num_parent = 0;
            for (int parent : built_neighbors[cur]) {
                int query_edge_idx = query_->GetEdgeIndex(parent, cur);
                for (int parent_cand : candidate_set_[parent]) {
                    for (int data_edge_idx : data_->GetIncidentEdges(parent_cand, cur_label)) {
                        int cand = data_->GetOppositePoint(data_edge_idx);
                        if (data_->GetDegree(cand) < query_->GetDegree(cur)) break;
                        if (data_->GetEdgeLabel(data_edge_idx) != query_->GetEdgeLabel(query_edge_idx)) continue;
                        if (num_visit_cs_[cand] < num_parent) continue;
                        if (data_->GetCoreNum(cand) < query_->GetCoreNum(cur)) continue;
                        if (num_visit_cs_[cand] == num_parent) {
                            num_visit_cs_[cand] += 1;
                            if (num_visit_cs_[cand] == 1) {
                                candidate_set_[cur].emplace_back(cand);
                                BitsetCS[cur][cand] = true;
                            }
                        }
                    }
                }
                num_parent++;
            }
            for (int j = 0; j < candidate_set_[cur].size(); j++) {
                int cand = candidate_set_[cur][j];
                BitsetCS[cur][cand] = false;
                if (num_visit_cs_[cand] == num_parent) {
                    BitsetCS[cur][cand] = true;
                }
                else {
                    candidate_set_[cur][j] = candidate_set_[cur].back();
                    candidate_set_[cur].pop_back();
                    j--;
                }
                num_visit_cs_[cand] = 0;
            }
            if (candidate_set_[cur].empty()) exit(2);
            for (int uc : query_->GetNeighbors(cur)) {
                if (candidate_set_[uc].empty()) {
                    built_neighbors[uc].push_back(cur);
                }
            }
        }

        int cs_edge = 0, cs_vertex = 0;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            for (int q_edge_idx : query_->GetAllIncidentEdges(i)) {
                int q_nxt = query_->GetOppositePoint(q_edge_idx);
                for (int j : candidate_set_[i]) {
                    for (int d_edge_idx : data_->GetIncidentEdges(j, query_->GetVertexLabel(q_nxt))) {
                        if (data_->GetEdgeLabel(d_edge_idx)!= query_->GetEdgeLabel(q_edge_idx)) continue;
                        int d_nxt = data_->GetOppositePoint(d_edge_idx);
                        if (BitsetCS[q_nxt][d_nxt]) {
                            BitsetEdgeCS[q_edge_idx][d_edge_idx] = true;
                            cs_edge++;
                        }
                    }
                }
            }
            cs_vertex += candidate_set_[i].size();
        }

        // for (int i = 0; i < query_->GetNumVertices(); i++) {
        //     fprintf(stderr, "Init CS %d: %lu\n", i, candidate_set_[i].size());
        // }
        return true;
    }

    bool CandidateSpace::InitRootCandidates(int root) {
        int root_label = query_->GetVertexLabel(root);
        for (int cand : data_->GetVerticesByLabel(root_label)) {
            if (data_->GetDegree(cand) < query_->GetDegree(root)) break;
            if (data_->GetCoreNum(cand) < query_->GetCoreNum(root)) continue;
            candidate_set_[root].emplace_back(cand);
            BitsetCS[root][cand] = true;
        }
        return !candidate_set_[root].empty();
    }

    bool CandidateSpace::RefineCS(){
        std::vector<int> local_stage(query_->GetNumVertices(), 0);
        std::vector<double> priority(query_->GetNumVertices(), 0.50);
        int queue_pop_count = 0;
        int maximum_queue_cnt = 5 * query_->GetNumEdges();
        int current_stage = 0;
        while (queue_pop_count < maximum_queue_cnt) {
            int cur = 0;
            double cur_priority = priority[cur];
            for (int i = 1; i < query_->GetNumVertices(); i++) {
                if (priority[i] > cur_priority) {
                    cur_priority = priority[i];
                    cur = i;
                }
                else if (priority[i] == cur_priority) {
                    if (GetCandidateSetSize(i) > GetCandidateSetSize(cur)) {
                        cur = i;
                    }
                }
            }
            if (cur_priority < 0.05) break;
            current_stage++;
            queue_pop_count+=query_->GetDegree(cur);
            int bef_cand_size = candidate_set_[cur].size();
            if (opt.neighborhood_filter == NEIGHBOR_SAFETY) {
                std::fill(neighbor_label_frequency.begin(), neighbor_label_frequency.end(), 0);
                memset(in_neighbor_cs, false, data_->GetNumVertices());
                PrepareNeighborSafety(cur);
            }
            for (int i = 0; i < candidate_set_[cur].size(); i++) {
                int cand = candidate_set_[cur][i];
                bool valid = true;
                if (valid and opt.structure_filter > NO_STRUCTURE_FILTER) valid = CheckSubStructures(cur, cand);
                if (valid) valid = NeighborFilter(cur, cand);
                // if (cur == 11 and cand == 56230) {
                //     fprintf(stderr, "Cur, Cand: %d, %d\n", cur, cand);
                //     fprintf(stderr, "NBR: %d\n", NeighborFilter(cur, cand));
                //     fprintf(stderr, "Struct: %d\n", CheckSubStructures(cur, cand));
                // }
                if (!valid) {
                    int removed = candidate_set_[cur][i];
                    for (int query_edge_idx : query_->GetAllIncidentEdges(cur)) {
                        int nxt = query_->GetOppositePoint(query_edge_idx);
                        int nxt_label = query_->GetVertexLabel(nxt);
                        for (int data_edge_idx : data_->GetIncidentEdges(removed, nxt_label)) {
                            int nxt_cand = data_->GetOppositePoint(data_edge_idx);
                            if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                            if (BitsetEdgeCS[query_edge_idx][data_edge_idx]) {
                                BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
                                BitsetEdgeCS[query_->GetOppositeEdge(query_edge_idx)][data_->GetOppositeEdge(data_edge_idx)] = false;
                            }
                        }
                    }
                    // if (cur == 11 and cand == 56230) {
                    //     fprintf(stderr, "Remove %d from CS[%d]! swap with %d\n",
                    //         candidate_set_[cur][i], cur, candidate_set_[cur].back());
                    // }
                    candidate_set_[cur][i] = candidate_set_[cur].back();
                    candidate_set_[cur].pop_back();
                    --i;
                    BitsetCS[cur][cand] = false;
                }
            }
            if (candidate_set_[cur].empty()) {
                exit(2);
            }
            int aft_cand_size = candidate_set_[cur].size();
            // fprintf(stderr, "Refine CS %d: %d -> %d\n", cur, bef_cand_size, aft_cand_size);

            if (aft_cand_size == bef_cand_size) {
                priority[cur] = 0;
                continue;
            }
            double out_prob = 1 - aft_cand_size * 1.0 / bef_cand_size;
            priority[cur] = 0;
            for (int nxt : query_->GetNeighbors(cur)) {
                priority[nxt] = 1 - (1 - out_prob) * (1 - priority[nxt]);
                local_stage[nxt] = current_stage;
            }
        }
        return true;
    }

    void CandidateSpace::PrepareNeighborSafety(int cur) {
        for (int q_neighbor : query_->GetNeighbors(cur)){
            neighbor_label_frequency[query_->GetVertexLabel(q_neighbor)]++;
            for (int d_neighbor : candidate_set_[q_neighbor]) {
                in_neighbor_cs[d_neighbor] = true;
            }
        }
    }

    bool CandidateSpace::CheckNeighborSafety(int cur, int cand) {
        for (int d_neighbor : data_->GetNeighbors(cand)) {
            if (in_neighbor_cs[d_neighbor]) {
                neighbor_label_frequency[data_->GetVertexLabel(d_neighbor)]--;
            }
        }
        bool valid = true;
        for (int l = 0; l < data_->GetNumLabels(); ++l) {
            if (neighbor_label_frequency[l] > 0) {
                valid = false;
                break;
            }
        }
        for (int d_neighbor : data_->GetNeighbors(cand)) {
            if (in_neighbor_cs[d_neighbor]) {
                neighbor_label_frequency[data_->GetVertexLabel(d_neighbor)]++;
            }
        }
        return valid;
    }

    inline bool CandidateSpace::EdgeCandidacy(int query_edge_id, int data_edge_id) {
        if (query_edge_id == -1 || data_edge_id == -1) {
            return false;
        }
        return BitsetEdgeCS[query_edge_id][data_edge_id];
    }

    inline bool CandidateSpace::TriangleSafety(int query_edge_id, int data_edge_id) {
        auto &query_triangles = query_->GetLocalTriangles(query_edge_id);
        if (query_triangles.empty()) return true;
        auto &candidate_triangles = data_->GetLocalTriangles(data_edge_id);
        if (query_triangles.size() > candidate_triangles.size()) return false;
        for (auto qtv: query_triangles) {
            bool found = std::any_of(candidate_triangles.begin(),
                                     candidate_triangles.end(),
                                     [&](auto tv) {
                                         return BitsetEdgeCS[get<1>(qtv)][get<1>(tv)] and  BitsetEdgeCS[get<2>(qtv)][get<2>(tv)];
                                     });
            if (!found) return false;
        }
        return true;
    };

    bool CandidateSpace::FourCycleSafety(int query_edge_id, int data_edge_id) {
        auto &query_cycles = query_->GetLocalFourCycles(query_edge_id);
        auto &data_cycles = data_->GetLocalFourCycles(data_edge_id);
        if (query_cycles.size() > data_cycles.size()) return false;
        for (int i = 0; i < query_cycles.size(); i++) {
            auto &q_info = query_cycles[i];
            for (int j = 0; j < data_cycles.size(); j++) {
                auto &d_info =  data_cycles[j];
                bool validity = true;
                validity &= BitsetEdgeCS[get<1>(q_info.edges)][get<1>(d_info.edges)];
                if (!validity) continue;
                validity &= BitsetEdgeCS[get<2>(q_info.edges)][get<2>(d_info.edges)];
                if (!validity) continue;
                validity &= BitsetEdgeCS[get<3>(q_info.edges)][get<3>(d_info.edges)];
                if (validity and get<0>(q_info.diags) != -1)
                    validity &= EdgeCandidacy(get<0>(q_info.diags),get<0>(d_info.diags));
                if (validity and get<1>(q_info.diags) != -1)
                    validity &= EdgeCandidacy(get<1>(q_info.diags), get<1>(d_info.diags));
                if (validity) {
                    goto nxt_cycle;
                }
            }
            return false;
            nxt_cycle:
            continue;
        }
        return true;
    };

    bool CandidateSpace::CheckSubStructures(int cur, int cand) {
        for (int query_edge_idx : query_->GetAllIncidentEdges(cur)) {
            int nxt = query_->GetOppositePoint(query_edge_idx);
            int nxt_label = query_->GetVertexLabel(nxt);
            bool found = false;
            for (int data_edge_idx : data_->GetIncidentEdges(cand, nxt_label)) {
                int nxt_cand = data_->GetOppositePoint(data_edge_idx);
                if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                if (!BitsetEdgeCS[query_edge_idx][data_edge_idx]) continue;
                if (!StructureFilter(query_edge_idx, data_edge_idx)) {
                    BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
                    BitsetEdgeCS[query_->GetOppositeEdge(query_edge_idx)][data_->GetOppositeEdge(data_edge_idx)] = false;
                    continue;
                }
                found = true;
            }
            if (!found) {
                return false;
            }
        }
        return true;
    };

    bool CandidateSpace::NeighborBipartiteSafety(int cur, int cand){
        if (query_->GetDegree(cur) == 1) {
            int uc = query_->GetNeighbors(cur)[0];
            int query_edge_index = query_->GetEdgeIndex(cur, uc);
            int label = query_->GetVertexLabel(uc);
            return std::any_of(data_->GetIncidentEdges(cand, label).begin(),
                               data_->GetIncidentEdges(cand, label).end(),
                               [&](int data_edge_index) {
                                   return BitsetEdgeCS[query_edge_index][data_edge_index];
                               });
        }
        BPSolver.Reset();
        int i = 0, j = 0;
        for (int &query_edge_index : query_->GetAllIncidentEdges(cur)) {
            j = 0;
            for (int &edge_id : data_->GetAllIncidentEdges(cand)) {
                if (BitsetEdgeCS[query_edge_index][edge_id]) {
                    BPSolver.AddEdge(i, j);
                }
                j++;
            }
            i++;
        }
        return BPSolver.Solve() == query_->GetDegree(cur);
    };

    bool CandidateSpace::EdgeBipartiteSafety(int cur, int cand) {
        // if (!(cur == 11 and cand == 56230)) {
        //     BPSolver.print_log = false;
        // }
        auto &query_edges = query_->GetAllIncidentEdges(cur);
        auto &data_edges = data_->GetAllIncidentEdges(cand);
        if (query_edges.size() == 1) {
            int q_edge_id = query_edges[0];
            for (int d_edge_id : data_edges) {
                if (BitsetEdgeCS[q_edge_id][d_edge_id])
                    return true;
            }
            return false;
        }
        std::vector<std::pair<int, int>> edge_pairs;
        int ii = 0, jj = 0;
        BPSolver.Reset();
        for (int query_edge_index : query_edges) {
            int uc = query_->GetOppositePoint(query_edge_index);
            jj = 0;
            for (int edge_id : data_edges) {
                int vc = data_->GetOppositePoint(edge_id);
                if (data_->GetDegree(vc) < query_->GetDegree(uc)) break;
                if (BitsetEdgeCS[query_edge_index][edge_id]) {
                    BPSolver.AddEdge(ii, jj);
                    edge_pairs.emplace_back(ii, jj);
                }
                jj++;
            }
            ii++;
        }
        // if (cur == 11 and cand == 56230) {
        //     printf("Checking %d, %d EBSafety\n", cur, cand);
        //     BPSolver.print_log = true;
        // }
        bool b = BPSolver.FindUnmatchableEdges(query_edges.size());
        if (!b) {
            // if (cur == 11 and cand == 56230) {
            //     printf("BP Matching Failed!\n");
            // }
            return false;
        }
        for (auto &[i, j] : edge_pairs) {
            if (!BPSolver.matchable[i][j]) {
                int left_unmatch = query_edges[i];
                int right_unmatch = data_edges[j];
                BitsetEdgeCS[left_unmatch][right_unmatch] = false;
                BitsetEdgeCS[query_->GetOppositeEdge(left_unmatch)][data_->GetOppositeEdge(right_unmatch)] = false;
            }
        }
        return true;
    }
}