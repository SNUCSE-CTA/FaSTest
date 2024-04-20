#pragma once
/**
* @brief Class for data graph in pattern matching related problems
* @type Assume undirected, connected, labeled graph
*/
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include "DataStructure/Graph.h"

namespace GraphLib { namespace SubgraphMatching {
    struct LabelStatistics {
        std::vector<double> vertex_label_probability, edge_label_probability;
        double vertex_label_entropy = 0.0, edge_label_entropy = 0.0;
    };
    class DataGraph : public GraphLib::Graph {
    protected:
        // array of vertices, grouped by label, ordered by decreasing order of degree
        std::vector<std::vector<int>> vertex_by_labels;
        // num_vertex_by_label_degree;
        std::unordered_map<int, int> transferred_label_map;
        LabelStatistics label_statistics;
    public:
        DataGraph(const Graph &g) : Graph(g) {};
        DataGraph(){};
        std::vector<int>& GetVerticesByLabel(int label) { return vertex_by_labels[label]; }
        inline int GetTransferredLabel(int l) {return transferred_label_map[l];}
        void Preprocess();
        void TransformLabel();
        void ComputeLabelStatistics();
        bool FourCycleEnumerated() {return !local_four_cycles.empty();}
    };



    void DataGraph::TransformLabel() {
        int cur_transferred_label = 0;
        for (int v = 0; v < GetNumVertices(); v++) {
            int l = vertex_label[v];
            if (transferred_label_map.find(l) == transferred_label_map.end()) {
                transferred_label_map[l] = cur_transferred_label;
                cur_transferred_label += 1;
            }
            vertex_label[v] = transferred_label_map[l];
            num_vertex_labels = std::max(num_vertex_labels, vertex_label[v]+1);
        }
    }


    void DataGraph::ComputeLabelStatistics() {
        label_statistics.vertex_label_probability.resize(GetNumLabels(), 1e-4);
        for (int i = 0; i < GetNumVertices(); i++) {
            label_statistics.vertex_label_probability[GetVertexLabel(i)] += 1.0;
        }
        for (int i = 0; i < GetNumLabels(); i++) {
            label_statistics.vertex_label_probability[i] /= (1.0 * GetNumVertices());
        }
        for (auto x : label_statistics.vertex_label_probability) {
            label_statistics.vertex_label_entropy -= x * log2(x);
        }
    }


    void DataGraph::Preprocess() {
        for (auto &it : adj_list) max_degree = std::max(max_degree, (int)it.size());
        TransformLabel();
        BuildIncidenceList();
        ComputeCoreNum();
        ComputeLabelStatistics();
        vertex_by_labels.resize(GetNumLabels());
        // num_vertex_by_label_degree.resize(GetNumLabels());
        for (int i = 0; i < GetNumVertices(); i++) {
            vertex_by_labels[GetVertexLabel(i)].push_back(i);
        }
        for (int i = 0; i < GetNumLabels(); i++) {
            if (vertex_by_labels[i].empty()) continue;
            std::stable_sort(vertex_by_labels[i].begin(), vertex_by_labels[i].end(), [this](int a, int b) {
                return GetDegree(a) > GetDegree(b);
            });
        }
        // std::cerr << "#Vertex = " << GetNumVertices() << " " << "#Edges = " << GetNumEdges() << std::endl;
        // std::cerr << "Degeneracy = " << degeneracy << std::endl;
        // std::cerr << "Ent(VertexLabel) = " << label_statistics.vertex_label_entropy << std::endl;
    }
} }
