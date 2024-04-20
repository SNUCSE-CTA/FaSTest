#pragma once
#include <unordered_map>
#include <map>
#include <iostream>
#include <fstream>
#include "Base/Base.h"
// #define HUGE_GRAPH
namespace GraphLib {
    class Graph {
    protected:
        // Adjacency List
        std::vector<std::vector<int>> adj_list;
        std::vector<int> core_num, vertex_color;
        std::vector<int> degeneracy_order;
        int num_vertex = 0, num_edge = 0, max_degree = 0, degeneracy = 0, num_color = 0;
        int num_vertex_labels = 0;

        /**
         * Basic data structures for graph
         * @attribute (vertex/edge)_label : array of labels
         * @attribute edge_list : list of edges as pair<int, int> form
         * @attribute incident_edges[v][l] : list of indices of incident edges from v and endpoint label l
         * @attribute all_incident_edges[v] : list of indices of all incident edges from v
         * @attribute edge_index_map : Queryable as map[{u, v}].
         */
        std::vector<int> vertex_label, edge_label, edge_to;
        std::vector<std::pair<int, int>> edge_list;
        std::vector<std::vector<int>> all_incident_edges;
#ifdef HUGE_GRAPH
        std::vector<std::map<int, std::vector<int>>> incident_edges;
#else
        std::vector<std::vector<std::vector<int>>> incident_edges;
#endif
        std::vector<std::unordered_map<int, int>> edge_index_map;


        /*
         * Enumeration of Small Cycles for Cyclic Substructure Filter
         * For each edge e, store local triangles and four-cycles
         */
        struct FourMotif {
            std::tuple<int, int, int, int> edges;
            std::tuple<int, int> diags;
            FourMotif(std::tuple<int, int, int, int> edges, std::tuple<int, int> diags) : edges(edges), diags(diags) {}
        };
        std::vector<std::vector<std::tuple<int, int, int>>> local_triangles;
        std::vector<std::vector<FourMotif>> local_four_cycles;
    public:
        Graph() {}
        ~Graph() {}
        Graph &operator=(const Graph &) = delete;

        std::vector<int>& GetNeighbors(int v) {
            return adj_list[v];
        }

        inline int GetDegree(int v) const {
            return adj_list[v].size();
        }

        inline int GetNumVertices() const {
            return num_vertex;
        }

        inline int GetNumEdges() const {
            return num_edge;
        }

        inline int GetMaxDegree() const {
            return max_degree;
        }

        void ComputeCoreNum();

        inline int GetCoreNum(int v) const {
            return core_num[v];
        }

        inline int GetDegeneracy() const {return degeneracy;}

        void AssignVertexColor();

        inline int GetNumColors() const {return num_color;}

        inline int GetVertexColor(int v) const {return vertex_color[v];}


        void LoadLabeledGraph(const std::string &filename);

        inline std::vector<int>& GetAllIncidentEdges(int v) {return all_incident_edges[v];}
        inline std::vector<int>& GetIncidentEdges(int v, int label) {return incident_edges[v][label];}
        inline int GetVertexLabel(int v) const {return vertex_label[v];}
        inline int GetEdgeLabel(int edge_id) const {return edge_label[edge_id];}
        inline int GetNumLabels() const {return num_vertex_labels;}
        inline int GetOppositeEdge(int edge_id) const {return edge_id^1;}
        inline int GetOppositePoint(int edge_id) const {return edge_to[edge_id];}
        inline int GetEdgeIndex(int u, int v) {
            auto it = edge_index_map[u].find(v);
            return (it == edge_index_map[u].end() ? -1 : it->second);
        }
        inline std::pair<int, int>& GetEdge(int edge_id) {
            return edge_list[edge_id];
        }

        inline std::vector<std::tuple<int, int, int>>& GetLocalTriangles(int edge_id) {return local_triangles[edge_id];}
        inline std::vector<FourMotif>& GetLocalFourCycles(int edge_id) {return local_four_cycles[edge_id];}

        void EnumerateLocalTriangles();
        void EnumerateLocalFourCycles();
        void ChibaNishizeki();

        /**
         * @brief Build the incidence list structure
         */
        void BuildIncidenceList();

        void LoadGraph(std::vector<int> &vertex_labels, std::vector<std::pair<int, int>> &edges,
                       std::vector<int> &edge_labels,   bool directed = false);

        void WriteToFile(string filename);
    };

    /**
     * @brief Compute the core number of each vertex
     * @date Oct 21, 2022
     */
    void Graph::ComputeCoreNum() {
        core_num.resize(num_vertex, 0);
        int *bin = new int[GetMaxDegree() + 1];
        int *pos = new int[GetNumVertices()];
        int *vert = new int[GetNumVertices()];

        std::fill(bin, bin + (GetMaxDegree() + 1), 0);

        for (int v = 0; v < GetNumVertices(); v++) {
            core_num[v] = adj_list[v].size();
            bin[core_num[v]] += 1;
        }

        int start = 0;
        int num;

        for (int d = 0; d <= GetMaxDegree(); d++) {
            num = bin[d];
            bin[d] = start;
            start += num;
        }

        for (int v = 0; v < GetNumVertices(); v++) {
            pos[v] = bin[core_num[v]];
            vert[pos[v]] = v;
            bin[core_num[v]] += 1;
        }

        for (int d = GetMaxDegree(); d--;)
            bin[d + 1] = bin[d];
        bin[0] = 0;

        for (int i = 0; i < GetNumVertices(); i++) {
            int v = vert[i];

            for (int u : GetNeighbors(v)) {
                if (core_num[u] > core_num[v]) {
                    int du = core_num[u];
                    int pu = pos[u];
                    int pw = bin[du];
                    int w = vert[pw];

                    if (u != w) {
                        pos[u] = pw;
                        pos[w] = pu;
                        vert[pu] = w;
                        vert[pw] = u;
                    }

                    bin[du]++;
                    core_num[u]--;
                }
            }
        }
        degeneracy_order.resize(GetNumVertices());
        for (int i = 0; i < GetNumVertices(); i++) {
            degeneracy_order[i] = vert[i];
        }
        std::reverse(degeneracy_order.begin(),degeneracy_order.end());

        degeneracy = 0;
        for (int i = 0; i < GetNumVertices(); i++) {
            degeneracy = std::max(core_num[i], degeneracy);
        }

        delete[] bin;
        delete[] pos;
        delete[] vert;
    }

    /**
     * @brief Greedy coloring of the graph, following the given initial order of vertices.
     * @date Sep 16, 2022
     */
    void Graph::AssignVertexColor() {
        vertex_color.resize(GetNumVertices(), -1);
        num_color = 0;
        bool *used = new bool[GetNumVertices()];
        for (int vertexID : degeneracy_order) {
            for (int neighbor : adj_list[vertexID]) {
                if (vertex_color[neighbor] == -1) continue;
                used[vertex_color[neighbor]] = true;
            }
            int c = 0; while (used[c]) c++;
            vertex_color[vertexID] = c;
            num_color = std::max(num_color, c+1);
            for (int neighbor : adj_list[vertexID]) {
                if (vertex_color[neighbor] == -1) continue;
                used[vertex_color[neighbor]] = false;
            }
        }
    }



    void Graph::LoadLabeledGraph(const std::string &filename) {
        std::ifstream fin(filename);
        int v, e;
        std::string ignore, type, line;
        fin >> ignore >> v >> e;
        num_vertex = v;
        // add edges in both directions
        num_edge = e * 2;
        adj_list.resize(num_vertex);
        vertex_label.resize(num_vertex);
        edge_label.resize(num_edge);
        int num_lines = 0;
        while (getline(fin, line)) {
            auto tok = parse(line, " ");
            type = tok[0];
            tok.pop_front();
            if (type[0] == 'v') {
                int id = std::stoi(tok.front());
                tok.pop_front();
                int l;
                if (tok.empty()) l = 0;
                else {
                    l = std::stoi(tok.front());
                    tok.pop_front();
                }
                vertex_label[id] = l;
            }
            else if (type[0] == 'e') {
                int v1, v2;
                v1 = std::stoi(tok.front()); tok.pop_front();
                v2 = std::stoi(tok.front()); tok.pop_front();
                adj_list[v1].push_back(v2);
                adj_list[v2].push_back(v1);
                edge_to.push_back(v2); edge_to.push_back(v1);
                edge_list.push_back({v1, v2});
                edge_list.push_back({v2, v1});
                int el = tok.empty() ? 0 : std::stoi(tok.front());
                edge_label[edge_list.size()-2] = edge_label[edge_list.size()-1] = el;
                max_degree = std::max(max_degree, (int)std::max(adj_list[v1].size(), adj_list[v2].size()));
            }
            num_lines++;
        }
    }

    void Graph::LoadGraph(std::vector<int> &vertex_labels,
                          std::vector<std::pair<int, int>> &edges,
                          std::vector<int> &edge_labels,
                          bool directed) {
        num_vertex = vertex_labels.size();
        num_edge = edges.size();
        if (!directed) num_edge *= 2;
        adj_list.resize(num_vertex);
        vertex_label.resize(num_vertex);
        edge_label.resize(num_edge);
        for (int i = 0; i < num_vertex; i++) vertex_label[i] = vertex_labels[i];
        for (int i = 0; i < edges.size(); i++) {
            auto &[v1, v2] = edges[i];
            int el = (edge_labels.size() > i) ? edge_labels[i] : 0;
            adj_list[v1].push_back(v2);
            edge_to.push_back(v2);
            edge_list.push_back({v1, v2});
            edge_label[edge_list.size()-1] = el;
//            fprintf(stderr, "Edge %d %d\n", v1, v2);
            if (!directed) {
                adj_list[v2].push_back(v1);
                edge_to.push_back(v1);
                edge_list.push_back({v2, v1});
                edge_label[edge_list.size()-2] = edge_label[edge_list.size()-1] = el;
            }
        }
    }


    void Graph::BuildIncidenceList() {
        all_incident_edges.resize(num_vertex);
        incident_edges.resize(num_vertex);
        edge_index_map.resize(num_vertex);
        for (int i = 0; i < GetNumVertices(); i++) {
#ifndef HUGE_GRAPH
            incident_edges[i].resize(GetNumLabels());
#endif
        }
        int edge_id = 0;
        for (auto& [u, v] : edge_list) {
            all_incident_edges[u].push_back(edge_id);
            incident_edges[u][GetVertexLabel(v)].push_back(edge_id);
            edge_index_map[u][v] = edge_id;
            edge_id++;
        }

        // sort edges by degree of endpoint
        for (int i = 0; i < GetNumVertices(); i++) {
#ifdef HUGE_GRAPH
            for (auto &[l, vec] : incident_edges[i]) {
                std::stable_sort(vec.begin(), vec.end(),[this](auto &a, auto &b) -> bool {
                    int opp_a = edge_list[a].second;
                    int opp_b = edge_list[b].second;
                    return adj_list[opp_a].size() > adj_list[opp_b].size();
                });
            }
#else
            for (auto &vec : incident_edges[i]) {
                std::stable_sort(vec.begin(), vec.end(),[this](auto &a, auto &b) -> bool {
                    int opp_a = edge_list[a].second;
                    int opp_b = edge_list[b].second;
                    return adj_list[opp_a].size() > adj_list[opp_b].size();
                });
            }
#endif
            std::stable_sort(all_incident_edges[i].begin(), all_incident_edges[i].end(), [this](auto &a, auto &b) -> bool {
                return adj_list[edge_list[a].second].size() > adj_list[edge_list[b].second].size();
            });
        }
    }

    void Graph::WriteToFile(std::string filename) {
        std::filesystem::path filepath = filename;
        std::filesystem::create_directories(filepath.parent_path());
        std::ofstream out(filename);
        out << "t " << GetNumVertices() << ' ' << GetNumEdges()/2 << '\n';
        for (int i = 0; i < GetNumVertices(); i++) {
            out << "v " << i << ' ' << GetVertexLabel(i) << ' ' << GetDegree(i) << '\n';
        }
        int idx = 0;
        for (auto &e : edge_list) {
            if (e.first < e.second) {
                out << "e " << e.first << ' ' << e.second << ' ' << GetEdgeLabel(idx) << '\n';
            }
            idx++;
        }
    }
}
