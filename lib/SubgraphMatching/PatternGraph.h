#pragma once
/**
* @brief Class for subgraph pattern
*/
#include "DataStructure/Graph.h"
#include "SubgraphMatching/DataGraph.h"

//#include "ortools/linear_solver/linear_solver.h"
//
//namespace OR = operations_research;
//using OR::MPSolver, OR::MPVariable, OR::MPConstraint;
//std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("PDLP"));
//const double infinity = solver->infinity();

namespace GraphLib::SubgraphMatching {
    class PatternGraph : public Graph {
    public:
        PatternGraph(){};
        PatternGraph(const Graph& g) : Graph(g) {};
        ~PatternGraph(){};

        PatternGraph &operator=(const PatternGraph &) = delete;
        PatternGraph(const PatternGraph &) = delete;

        void ProcessPattern(DataGraph &data);
        std::vector<std::vector<int>> adj_idx;
        inline const int GetAdjIdx(int u, int uc) {return adj_idx[u][uc];};

//        void FindFractionalEdgeCover(std::vector<double> &weights);
//
//        std::vector<double> fractional_edge_cover;
    };

    void PatternGraph::ProcessPattern(DataGraph &data) {
        // transfer label & construct adj list and label frequency
        max_degree = 0;
        for (int v = 0; v < GetNumVertices(); v++) {
            int l = data.GetTransferredLabel(GetVertexLabel(v));
            vertex_label[v] = l;
            max_degree = std::max(max_degree, (int)(adj_list[v].size()));
        }
        num_vertex_labels = data.GetNumLabels();
        BuildIncidenceList();
        ComputeCoreNum();
        adj_idx.resize(GetNumVertices(), std::vector<int>(GetNumVertices(), -1));
        for (int u = 0; u < GetNumVertices(); u++) {
            for (int i = 0; i < adj_list[u].size(); i++) {
                adj_idx[u][adj_list[u][i]] = i;
            }
        }
    }

    /*void PatternGraph::FindFractionalEdgeCover(std::vector<double> &weights) {
        std::vector<OR::MPVariable*> edgevariables(GetNumEdges()/2);
        std::vector<OR::MPConstraint*> vertexconstraints(GetNumVertices());
        OR::MPObjective* objective = solver->MutableObjective();
        std::vector<int> X[50];
        for (int i = 0; i < GetNumVertices(); i++) {
            vertexconstraints[i] = solver->MakeRowConstraint(1.0, infinity);
        }
        for (int i = 0; i < GetNumEdges(); i+=2) {
            edgevariables[i/2] = solver->MakeNumVar(0.0, 1.0, "e" + std::to_string(i/2));
            auto &[u, v] = edge_list[i];
            vertexconstraints[u]->SetCoefficient(edgevariables[i/2], 1.0);
            vertexconstraints[v]->SetCoefficient(edgevariables[i/2], 1.0);
            X[u].push_back(i/2);
            X[v].push_back(i/2);
//            std::cout << "weight " << i/2 << " : " << weights[i] << std::endl;
            objective->SetCoefficient(edgevariables[i/2], weights[i]);
        }
//        for (int i = 0; i < GetNumVertices(); i++) {
//            printf("Constraint [%d]: ", i);
//            for (auto &it : X[i]) {
//                printf(" %d ", it);
//            }
//            printf("\n");
//        }
        objective->SetMinimization();
        const MPSolver::ResultStatus result_status = solver->Solve();
        // Check that the problem has an optimal solution.
//        if (result_status != MPSolver::OPTIMAL) {
//            LOG(FATAL) << "The problem does not have an optimal solution!";
//        }

        fractional_edge_cover.resize(GetNumEdges()/2);
//        LOG(INFO) << "Solution:";
//        LOG(INFO) << "Optimal objective value = " << objective->Value();
        for (int i = 0; i < GetNumEdges(); i+=2) {
//            LOG(INFO) << edgevariables[i/2]->name() << " = " << edgevariables[i/2]->solution_value();
            fractional_edge_cover[i/2] = edgevariables[i/2]->solution_value();
        }
    }*/
}
