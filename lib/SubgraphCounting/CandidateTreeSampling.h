#pragma once
#include <boost/math/distributions.hpp>
#include "SubgraphMatching/CandidateSpace.h"

static std::random_device rd;
static std::mt19937 gen(rd());
inline int sample_from_distribution(std::discrete_distribution<int> &weighted_distr) {
    return weighted_distr(gen);
}

namespace GraphLib {
    using SubgraphMatching::DataGraph, SubgraphMatching::PatternGraph;
    using SubgraphMatching::CandidateSpace;
    using std::vector;
    namespace CardinalityEstimation {
        struct QueryTree {
            PatternGraph *query_;
            vector<vector<int>> tree_adj_list, tree_children;
            vector<int> tree_sequence;
            vector<int> parent, child_index;
            int root;

            void Initialize(PatternGraph *query, int root_idx) {
                query_ = query;
                tree_adj_list.clear();
                tree_children.clear();
                tree_sequence.clear();
                parent.clear();
                child_index.clear();
                root = root_idx;

                tree_adj_list.resize(query_->GetNumVertices());
                tree_children.resize(query_->GetNumVertices());
                tree_sequence.resize(query_->GetNumVertices(), -1);
                parent.resize(query_->GetNumVertices(), -1);
                child_index.resize(query_->GetNumVertices(), -1);
            }

            void AddEdge(int u, int v) {
                tree_adj_list[u].push_back(v);
                tree_adj_list[v].push_back(u);
            }

            void BuildTree() {
                bool *visit = new bool[query_->GetNumVertices()];
                std::fill(visit, visit + query_->GetNumVertices(), false);
                std::queue<int> q;
                int id = 0;
                parent[root] = -1;
                tree_sequence[id++] = root;
                q.push(root);
                visit[root] = true;
                while (!q.empty()) {
                    int v = q.front();
                    q.pop();
                    for (int c : tree_adj_list[v]) {
                        if (!visit[c]) {
                            q.push(c);
                            visit[c] = true;
                            child_index[c] = tree_children[v].size();
                            tree_children[v].push_back(c);
                            parent[c] = v;
                            tree_sequence[id++] = c;
                        }
                    }
                }
                delete[] visit;
            }

            vector<int>& GetChildren(int v) {return tree_children[v];}
            int GetParent(int v) {return parent[v];}
            int GetChildIndex(int v) {return child_index[v];}
            int GetKthVertex(int k) {return tree_sequence[k];}
        };
        class CandidateTreeSampler {
        protected:
            dict info;
            DataGraph *data_;
            PatternGraph *query_;
            CandidateSpace *CS;
            CardEstOption opt;

            double **num_trees_, total_trees_;
            bool *seen;
            // vector<vector<vector<vector<int>>>> sample_candidates_;
            vector<vector<vector<vector<double>>>> sample_candidate_weights_;
            vector<vector<vector<std::discrete_distribution<int>>>> sample_dist_;
            vector<int> root_candidates_, sample;
            vector<double> root_weights_;
            std::discrete_distribution<int> sample_root_dist_;

            QueryTree Tq;
        public:
            dict GetInfo() {return info;}
            CandidateTreeSampler(DataGraph *data, CardEstOption option) {
                data_ = data;
                opt = option;
                num_trees_ = new double *[opt.MAX_QUERY_VERTEX];
                for (int i = 0; i < opt.MAX_QUERY_VERTEX; i++) {
                    num_trees_[i] = new double[data_->GetNumVertices()];
                }
                seen = new bool[data_->GetNumVertices()]();
            }


            void Preprocess(PatternGraph *query, CandidateSpace *cs) {
                info.clear();
                Timer timer; timer.Start();
                query_ = query;
                CS = cs;
                sample_dist_.clear();
                // sample_candidates_.clear();
                sample_candidate_weights_.clear();
                root_candidates_.clear();

                sample.resize(query_->GetNumVertices(), -1);
                std::memset(seen, 0, sizeof(bool) * data_->GetNumVertices());
                BuildSpanningTree();
                CountCandidateTrees();
                timer.Stop();
                info["TreeCountTime"] = timer.GetTime();
                info["#CandTree"] = total_trees_;
            };

            void BuildSpanningTree() {
                if (opt.treegen_strategy == CardinalityEstimation::TREEGEN_RANDOM) {
                    int root_node = rand() % query_->GetNumVertices();
                    Tq.Initialize(query_, root_node);
                    int num_discovered = 1; 
                    std::vector<int> is_discovered(query_->GetNumVertices(), false);
                    is_discovered[root_node] = true;
                    int cur = root_node;
                    while (num_discovered < query_->GetNumVertices()) {
                        int d = query_->GetDegree(cur);
                        int di = rand() % d;
                        int v = query_->GetNeighbors(cur)[di];
                        if (!is_discovered[v]) {
                            Tq.AddEdge(cur, v);
                            is_discovered[v] = true;
                            num_discovered++;
                        }
                        cur = v; 
                    }
                    Tq.BuildTree();
                    return; 
                }
                int root_idx = 0, num_root_cands = CS->GetCandidateSetSize(0);
                std::vector<std::pair<double, std::pair<int, int>>> edges;
                for (int i = 0; i < query_->GetNumVertices(); i++) {
                    for (int q_neighbor : query_->GetNeighbors(i)) {
                        if (i > q_neighbor) continue;
                        double density = 0.0;
                        for (int cand_idx = 0; cand_idx < CS->GetCandidateSetSize(i); cand_idx++) {
                            int num_cs_neighbor = CS->GetCandidateNeighbors(i, cand_idx, q_neighbor).size();
                            density += num_cs_neighbor;
                        }
                        if (opt.treegen_strategy == CardinalityEstimation::TREEGEN_DENSITY_MST) {
                            density /= ((CS->GetCandidateSetSize(i) * 1.0) * (CS->GetCandidateSetSize(q_neighbor) * 1.0));
                        }
                        if (density > 0) {
                            edges.push_back({density * 1.0, {i, q_neighbor}});
                        }
                    }
                    if (CS->GetCandidateSetSize(i) < num_root_cands) {
                        num_root_cands = CS->GetCandidateSetSize(i);
                        root_idx = i;
                    }
                }

                std::sort(edges.begin(), edges.end());
                Tq.Initialize(query_, root_idx);
                int num_tree_edges = 0;
                int qV = query_->GetNumVertices();
                std::vector<int> deg(qV, 0);
                UnionFind uf(qV);
                while (num_tree_edges + 1 < qV) {
                    double minw = 1e9;
                    std::pair<int, int> me;
                    for (auto &[w, e] : edges) {
                        auto [u, v] = e;
                        if (uf.find(u) == uf.find(v)) continue;
                        if (minw > w + deg[u] * 1e-7) {
                            minw = w + deg[u] * 1e-7;
                            me = e;
                        }
                    }
                    uf.unite(me.first, me.second);
                    deg[me.first]++;
                    deg[me.second]++;
                    Tq.AddEdge(me.first, me.second);
                    num_tree_edges++;
                }
                Tq.BuildTree();
            };


            void CountCandidateTrees() {
                for (int i = 0; i < query_->GetNumVertices(); i++) {
                    memset(num_trees_[i], 0, sizeof(double) * CS->GetCandidateSetSize(i));
                }
                int cnt = 0;
                sample_candidate_weights_.resize(query_->GetNumVertices());
                sample_dist_.resize(query_->GetNumVertices());
                for (int i = 0; i < query_->GetNumVertices(); i++) {
                    int u = Tq.GetKthVertex(query_->GetNumVertices() - i - 1);
                    int num_cands = CS->GetCandidateSetSize(u);
                    auto children = Tq.GetChildren(u);
                    int num_children = children.size();
                    sample_candidate_weights_[u].resize(num_cands);
                    sample_dist_[u].resize(num_cands);

                    std::vector<double> tmp_num_child(num_children);
                    for (int cs_idx = 0; cs_idx < num_cands; cs_idx++) {
                        sample_candidate_weights_[u][cs_idx].resize(num_children);
                        sample_dist_[u][cs_idx].resize(num_children);

                        double num_ = 1.0;
                        std::fill(tmp_num_child.begin(), tmp_num_child.end(), 0.0);
                        for (int uc_idx = 0; uc_idx < num_children; uc_idx++) {
                            int uc = children[uc_idx];
                            auto candidate_neighbors = CS->GetCandidateNeighbors(u, cs_idx, uc);
                            sample_candidate_weights_[u][cs_idx][uc_idx].resize(candidate_neighbors.size());
                            for (int j = 0; j < candidate_neighbors.size(); j++) {
                                int vc_idx = candidate_neighbors[j];
                                tmp_num_child[uc_idx] += num_trees_[uc][vc_idx];
                                sample_candidate_weights_[u][cs_idx][uc_idx][j] = num_trees_[uc][vc_idx];
                            }
                        }
                        for (int j = 0; j < num_children; j++) {
                            num_ *= tmp_num_child[j];
                            sample_dist_[u][cs_idx][j] = std::discrete_distribution<int>(
                                    sample_candidate_weights_[u][cs_idx][j].begin(),
                                    sample_candidate_weights_[u][cs_idx][j].end());
                            cnt++;
                        }
                        num_trees_[u][cs_idx] = num_;
                    }
                }

                total_trees_ = 0.0;
                root_candidates_.clear();
                root_weights_.clear();
                int root = Tq.root;
                int root_candidate_size = CS->GetCandidateSetSize(root);
                for (int root_candidate_idx = 0; root_candidate_idx < root_candidate_size; ++root_candidate_idx) {
                    total_trees_ += num_trees_[root][root_candidate_idx];
                    if (num_trees_[root][root_candidate_idx] > 0) {
                        root_candidates_.emplace_back(root_candidate_idx);
                        root_weights_.emplace_back(num_trees_[root][root_candidate_idx]);
                    }
                }
                sample_root_dist_ = std::discrete_distribution<int>(root_weights_.begin(), root_weights_.end());
            };

            bool GetSampleTree() {
                std::fill(sample.begin(), sample.end(), -1);
                bool valid = true;
                sample[Tq.root] = root_candidates_[sample_from_distribution(sample_root_dist_)];
                seen[CS->GetCandidate(Tq.root, sample[Tq.root])] = true;
                for (int i = 0; i < query_->GetNumVertices(); ++i) {
                    int u = Tq.GetKthVertex(i);
                    int v_idx = sample[u];
                    auto &children = Tq.GetChildren(u);
                    for (int uc_idx = 0; uc_idx < children.size(); ++uc_idx) {
                        int uc = children[uc_idx];
                        int vc_idx = sample_from_distribution(sample_dist_[u][v_idx][uc_idx]);
                        sample[uc] = CS->GetCandidateNeighbor(u, v_idx, uc, vc_idx);
                        int cand = CS->GetCandidate(uc, sample[uc]);
                        if (seen[cand]) {
                            valid = false;
                            goto INJECTIVE_VIOLATED;
                        }
                        seen[cand] = true;
                    }
                }
                INJECTIVE_VIOLATED:
                for (int i = 0; i < query_->GetNumVertices(); i++) {
                    if (sample[i] >= 0) {
                        int cand = CS->GetCandidate(i, sample[i]);
                        seen[cand] = false;
                    }
                }
                if (!valid) return false;
                for (int i = 0; i < query_->GetNumVertices(); i++) {
                    if (sample[i] == -1) return false;
                    sample[i] = CS->GetCandidate(i, sample[i]);
                }
                for (int i = 0; i < query_->GetNumVertices(); i++) {
                    for (int qe : query_->GetAllIncidentEdges(i)) {
                        int j = query_->GetOppositePoint(qe);
                        int de = data_->GetEdgeIndex(sample[i], sample[j]);
                        if (de == -1)
                            return false;
                        if (data_->GetEdgeLabel(de) != query_->GetEdgeLabel(qe))
                            return false;
                    }
                }
                return true;
            }

            // (Estimate, #Success)
            std::pair<double, int> Estimate() {
                Timer timer; timer.Start();
                int success = 0, trials = 0;
                while (++trials) {
                    auto result = GetSampleTree();
                    if (result) success++;
                    double rhohat = (success * 1.0 / trials);
                    if (trials == 50000 and success <= 10) {
                        timer.Stop();
                        info["#TreeTrials"] = trials;
                        info["#TreeSuccess"] = success;
                        info["TreeSampleTime"] = timer.GetTime();
                        return {-1, success};
                    }
                    if (trials >= 1000 and trials % 100 == 0) {
                        long double wplus = boost::math::binomial_distribution<>::find_upper_bound_on_p(trials, success, 0.05/2);
                        long double wminus = boost::math::binomial_distribution<>::find_lower_bound_on_p(trials, success, 0.05/2);
                        if (rhohat * 0.8 < wminus && wplus < rhohat * 1.25) {
                            timer.Stop();
                            break;
                        }
                    }
                }
                auto est = std::make_pair((success * 1.0 / (trials * 1.0)) * total_trees_, success);
                info["#TreeTrials"] = trials;
                info["#TreeSuccess"] = success;
                info["TreeSampleTime"] = timer.GetTime();
                return est;
            }
        };
    }
}