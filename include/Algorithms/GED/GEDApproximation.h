//
// Created by florian on 28.08.25.
//

#ifndef TESTGRAPHLIB_GEDAPPROXIMATION_H
#define TESTGRAPHLIB_GEDAPPROXIMATION_H
#include <iosfwd>

#include "GraphDataStructures/GraphBase.h"

#endif //TESTGRAPHLIB_GEDAPPROXIMATION_H



// enum of the type of approximation
enum class GEDApproximationType {
  ASTAR,
};

inline std::ostream& operator<<(std::ostream& os, const GEDApproximationType& type) {
  switch (type) {
    case GEDApproximationType::ASTAR:
      return os << "A*";
  }
  return os;
}

enum class EditType {
  INSERT,
  DELETE,
  REPLACE,
  RELABEL,
};

inline std::ostream& operator<<(std::ostream& os, const EditType& type) {
  switch (type) {
    case EditType::INSERT:
      return os << "INSERT";
    case EditType::DELETE:
      return os << "DELETE";
    case EditType::REPLACE:
      return os << "REPLACE";
    case EditType::RELABEL:
      return os << "RELABEL";
    default:
      return os;
  }
}


enum class OperationObject {
  NODE,
  EDGE,
};

inline std::ostream& operator<<(std::ostream& os, const OperationObject& operationObject) {
  switch (operationObject) {
    case OperationObject::NODE:
      return os << "NODE";
    case OperationObject::EDGE:
      return os << "EDGE";
    default:
      return os;
  }
}


struct EditOperation {
  OperationObject operationObject;
  EditType type;
  NodeId node;
  EDGE edge;
};

inline std::ostream& operator<<(std::ostream& os, const EditOperation& operation) {
  switch (operation.operationObject) {
    case OperationObject::NODE:
       os << "NODE " << operation.node;
      break;
    case OperationObject::EDGE:
       os << "EDGE " << operation.edge.first << "--" << operation.edge.second;
      break;
  }
  switch (operation.type) {
    case EditType::INSERT:
      os << "INSERT" << std::endl;
    case EditType::DELETE:
      os << "DELETE" << std::endl;
    case EditType::REPLACE:
      os << "REPLACE" << std::endl;
    case EditType::RELABEL:
      os << "RELABEL" << std::endl;
  }
  return os;
}


struct GEDApproximationParameters {

  struct Costs {
    double node_sub = 1.0;  // cost if node labels differ
    double node_del = 1.0;
    double node_ins = 1.0;
    double edge_del = 1.0;
    double edge_ins = 1.0;
  };

  // Input Parameters
  GEDApproximationType type = GEDApproximationType::ASTAR;
  Costs costs;
  int seed = 0;
  int max_expansions = 200000;
  int expansions = 0;
  std::vector<int> mapping;  // size n1, -1 means deleted; unmapped G2 are inserted
  bool complete = false;


  // Output Parameters
  double approximated_distance = std::numeric_limits<double>::infinity();
  GraphData<GraphStruct> edit_path_graphs;
  std::vector<EditOperation> edit_operations;
  double runtime = 0.0;



};




class GEDApproximation {
public:
  GEDApproximation(GraphStruct& source, GraphStruct& target) : _source(source), _target(target), _parameters() {
  };
  void Run(GEDApproximationParameters& parameters);


  void AStar(GEDApproximationParameters& parameters);

private:
  GraphStruct& _source;
  GraphStruct& _target;
  GEDApproximationParameters _parameters;

  struct State {
    // map12[i] in [0..n2-1] => node i (G1) mapped to node map12[i] (G2)
    // map12[i] == -1 => node i deleted (ε)
    std::vector<int> map12;
    std::vector<char> used2;     // size n2: whether a G2 node already used
    int next_i = 0;         // next G1 node index to assign
    double g = 0.0;         // cost so far
    double h = 0.0;         // heuristic
    // priority is f = g + h
  };

  struct PQCmp {
    bool operator()(const State& a, const State& b) const {
      double fa = a.g + a.h, fb = b.g + b.h;
      if (fa != fb) return fa > fb; // min-heap
      // tie-break: deeper first to reduce branching
      if (a.next_i != b.next_i) return a.next_i < b.next_i;
      return a.g > b.g;
    }
  };


  // Helper functions
  std::pair<double, std::vector<int>> hungarian_min(const std::vector<std::vector<double>>& matrix);
  double heuristic_LSAP(const std::vector<int>& map12, const std::vector<char>& used2, int next_i);
  std::pair<double, std::vector<int>> initial_upper_bound();
  double totalCost_given_full_phi(const std::vector<int>& phi) const;
  double nodeSubCost(int i, int j) const;
  double edgeMismatchCost(const bool e1, bool e2) const;
  double heuristic_for_state(const State& s);
  double incEdgeCost_assign_i_to_j(const std::vector<int>& map12, int i, int j) const;
};

inline void GEDApproximation::Run(GEDApproximationParameters &parameters) {
  parameters.edit_path_graphs.add(_source);
  switch (parameters.type) {
    case GEDApproximationType::ASTAR:
      AStar(parameters);
      break;
    default:
      break;
  }

}

inline void GEDApproximation::AStar(GEDApproximationParameters &parameters) {
  const int n1 = _source.nodes();
  const int n2 = _target.nodes();

  // Initial UB from node-only LSAP
  auto [best_cost, best_map] = initial_upper_bound();

  std::priority_queue<State, std::vector<State>, PQCmp> pq;

  State start;
  start.map12.assign(n1, -2); // -2 = unassigned; -1 = delete; >=0 = mapped
  start.used2.assign(n2, false);
  start.next_i = 0;
  start.g = 0.0;
  start.h = heuristic_LSAP(start.map12, start.used2, start.next_i);
  pq.push(start);

  // simple prefix DP to prune dominated states
  std::unordered_map<std::string, double> best_g_for_prefix;
  best_g_for_prefix.reserve(1<<12);

  auto prefix_key = [&](const std::vector<int>& m, int upto)->std::string{
      std::string key; key.reserve(upto*3 + 2);
      for (int i = 0; i < upto; ++i) {
          int v = m[i];
          if (v == -2) key += "u,"; // unreachable here
          else if (v == -1) key += "d,";
          else { key += '0' + (v%10); key += ','; } // compact-ish (fine for pruning)
      }
      return key;
  };

  _parameters.approximated_distance = best_cost;
  _parameters.mapping = best_map;
  _parameters.complete = false;
  _parameters.expansions = 0;

  while (!pq.empty()) {
      if (_parameters.expansions >= _parameters.max_expansions) break;

      State s = pq.top(); pq.pop();
      double f = s.g + s.h;
      if (f >= _parameters.approximated_distance - 1e-12) continue; // pruned by best UB

      // Goal check (all G1 nodes assigned)
      if (s.next_i >= n1) {
          // Convert to exact total and maybe update best
          std::vector<int> phi = s.map12;
          for (int i = 0; i < n1; ++i) if (phi[i] == -2) phi[i] = -1; // sanity
          double exact = totalCost_given_full_phi(phi);
          if (exact + 1e-12 < _parameters.approximated_distance) {
              _parameters.approximated_distance = exact;
              _parameters.mapping = phi;
              _parameters.complete = true; // we reached a concrete full mapping
          }
          continue;
      }

      // prefix pruning
      std::string key = prefix_key(s.map12, s.next_i);
      auto it = best_g_for_prefix.find(key);
      if (it != best_g_for_prefix.end() && s.g >= it->second - 1e-12) {
          continue;
      }
      best_g_for_prefix[key] = s.g;

      int i = s.next_i;

      // Branch: try mapping i -> each free j in G2
      for (int j = 0; j < n2; ++j) if (!s.used2[j]) {
          State t = s;
          t.map12[i] = j;
          t.used2[j] = true;
          t.next_i = i + 1;
          double c_node = nodeSubCost(i, j);
          double c_edge = incEdgeCost_assign_i_to_j(s.map12, i, j);
          t.g = s.g + c_node + c_edge;
          t.h = heuristic_for_state(t);
          if (t.g + t.h < _parameters.approximated_distance - 1e-12) {
              pq.push(std::move(t));
          }
      }

      // Branch: i -> ε (delete)
      {
          State t = s;
          t.map12[i] = -1;
          t.next_i = i + 1;
          double c_node = _parameters.costs.node_del;
          double c_edge = incEdgeCost_assign_i_to_j(s.map12, i, -1);
          t.g = s.g + c_node + c_edge;
          t.h = heuristic_for_state(t);
          if (t.g + t.h < _parameters.approximated_distance - 1e-12) {
              pq.push(std::move(t));
          }
      }

      ++_parameters.expansions;
  }
}

// Returns (minCost, assignment) where assignment[i] = j is the chosen column for row i.
inline std::pair<double, std::vector<int>> GEDApproximation::hungarian_min(const std::vector<std::vector<double>>& matrix) {
  // Implementation adapted from classic O(n^3) algorithm with potentials.
  // Assumes finite costs (>=0). Works with doubles.
  const int n = static_cast<int>(matrix.size());
  std::vector<double> u(n + 1, 0), v(n + 1, 0); // potentials
  std::vector<int> p(n + 1, 0), way(n + 1, 0);

  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    std::vector<double> min_value(n + 1, std::numeric_limits<double>::infinity());
    std::vector<char> used(n + 1, false);
    int j0 = 0;
    do {
      used[j0] = true;
      int i0 = p[j0], j1 = 0;
      double delta = std::numeric_limits<double>::infinity();
      for (int j = 1; j <= n; ++j) if (!used[j]) {
        double cur = matrix[i0-1][j-1] - u[i0] - v[j];
        if (cur < min_value[j]) { min_value[j] = cur; way[j] = j0; }
        if (min_value[j] < delta) { delta = min_value[j]; j1 = j; }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
        else          { min_value[j] -= delta; }
      }
      j0 = j1;
    } while (p[j0] != 0);
    do {
      int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }
  std::vector<int> assignment(n, -1);
  for (int j = 1; j <= n; ++j) if (p[j] != 0) assignment[p[j]-1] = j-1;
  double value = -v[0];
  return {value, assignment};
}


// Compute initial upper bound via a node-only LSAP, then score it exactly (including edges).
inline std::pair<double, std::vector<int>> GEDApproximation::initial_upper_bound() {
  const int n1 = _source.nodes();
  const int n2 = _target.nodes();
  const int r1 = n1;
  const int r2 = n2;
  const int k = std::max(r1, r2);
  std::vector<std::vector<double>> M(k, std::vector<double>(k, 0.0));
  // Rows: G1 nodes + ε, Cols: G2 nodes + ε
  for (int ri = 0; ri < k; ++ri) {
    bool row_is_node = (ri < n1);
    int i = row_is_node ? ri : -1;
    for (int cj = 0; cj < k; ++cj) {
      bool col_is_node = (cj < n2);
      int j = col_is_node ? cj : -1;
      double c = 0.0;
      if (row_is_node && col_is_node) {
        c = this->nodeSubCost(i, j);
      } else if (row_is_node && !col_is_node) {
        c = _parameters.costs.node_del;
      } else if (!row_is_node && col_is_node) {
        c = _parameters.costs.node_ins;
      } else c = 0.0;
      M[ri][cj] = c;
    }
  }
  auto [_, assign] = hungarian_min(M);

  std::vector<int> phi(n1, -1);
  for (int ri = 0; ri < k; ++ri) {
    int cj = assign[ri];
    bool row_is_node = (ri < n1);
    bool col_is_node = (cj < n2);
    if (row_is_node && col_is_node) phi[ri] = cj;
  }
  double ub = totalCost_given_full_phi(phi);
  return {ub, phi};
}

// Exact total cost for a complete node mapping phi : V(G1)->V(G2)∪{ε}.
// Unmapped G2 nodes are treated as insertions; includes *all* node and edge costs.
inline double GEDApproximation::totalCost_given_full_phi(const std::vector<int>& phi) const {
  const int n1 = _source.nodes();
  const int n2 = _target.nodes();
  std::vector<char> used2(n2, false);
  double cost = 0.0;

  // Node costs and record used2
  for (int i = 0; i < n1; ++i) {
    int j = phi[i];
    if (j >= 0) {
      used2[j] = true;
      cost += nodeSubCost(i, j);
    } else { // deletion
      cost += _parameters.costs.node_del;
    }
  }
  // Insert remaining G2 nodes
  std::vector<int> inserted;
  for (int j = 0; j < n2; ++j) if (!used2[j]) {
    cost += _parameters.costs.node_ins;
    inserted.push_back(j);
  }

  // Edge costs for pairs of G1 nodes
  for (int i = 0; i < n1; ++i) for (int p = i+1; p < n1; ++p) {
    bool e1 = _source.edge(i, p);
    int j = phi[i], r = phi[p];
    bool e2 = false;
    if (j >= 0 && r >= 0) e2 = _target.edge(j, r);
    // If at least one endpoint missing in G2 (deleted), edge in G2 is absent.
    cost += edgeMismatchCost(e1, e2);
  }

  // Edge insertions in G2 involving at least one inserted node.
  // (Pairs where both endpoints are images of some G1 nodes were handled above.)
  // Inserted ↔ mapped, and inserted ↔ inserted
  std::vector<int> mapped2;
  mapped2.reserve(n2 - (int)inserted.size());
  for (int j = 0; j < n2; ++j) if (used2[j]) mapped2.push_back(j);

  // inserted ↔ mapped
  for (int a : inserted) for (int b : mapped2) {
    if (_target.edge(a, b)) cost += _parameters.costs.edge_ins;
  }
  // inserted ↔ inserted
  for (int x = 0; x < (int)inserted.size(); ++x)
    for (int y = x+1; y < (int)inserted.size(); ++y)
      if (_target.edge(inserted[x], inserted[y])) cost += _parameters.costs.edge_ins;

  return cost;
}

inline double GEDApproximation::nodeSubCost(int i, int j) const {
  // check wheter both graphs are labeled
  if (_source.labelType == _target.labelType && _source.labelType != LABEL_TYPE::UNLABELED) {
    return (_source.label(i) == _target.label(j)) ? 0.0 : _parameters.costs.node_sub;
  }
  return 0.0;
}

inline double GEDApproximation::edgeMismatchCost(const bool e1, const bool e2) const {
  if (e1 && !e2) return _parameters.costs.edge_del;
  if (!e1 && e2) return _parameters.costs.edge_ins;
  return 0.0; // both same (0/0 or 1/1 with unlabeled edges)
}

// Build an LSAP (Hungarian) lower bound for a partial mapping state.
// Rows = remaining G1 nodes + ε-rows; Cols = remaining G2 nodes + ε-cols.
// Costs include (i) node sub/del/ins, plus (ii) *only* edge mismatches to already-assigned nodes.
// This is admissible; it ignores edges between as-yet-unassigned nodes.
inline double GEDApproximation::heuristic_LSAP(const std::vector<int>& map12, const std::vector<char>& used2, const int next_i) {
  const int n1 = _source.nodes();
    const int n2 = _target.nodes();

    // Gather remaining indices
    std::vector<int> rem1;
    for (int i = next_i; i < n1; ++i) rem1.push_back(i); // we assign in order
    std::vector<int> rem2;
    for (int j = 0; j < n2; ++j) if (!used2[j]) rem2.push_back(j);

    int r1 = (int)rem1.size();
    int r2 = (int)rem2.size();

    // Boundary cases where we can compute exact residual (still admissible).
    if (r1 == 0 && r2 == 0) return 0.0;

    if (r1 == 0) {
        // Only insertions left: exact node insert cost + edges to assigned + edges among remaining
        double h = 0.0;
        // to assigned
        for (int jj : rem2) {
            h += _parameters.costs.node_ins;
            for (int p = 0; p < n1; ++p) {
                int q = map12[p];
                if (q >= 0 && _target.edge(jj, q)) h += _parameters.costs.edge_ins;
            }
        }
        // among remaining new nodes
        for (int x = 0; x < (int)rem2.size(); ++x)
            for (int y = x+1; y < (int)rem2.size(); ++y)
                if (_target.edge(rem2[x], rem2[y])) h += _parameters.costs.edge_ins;
        return h;
    }

    if (r2 == 0) {
        // Only deletions left: exact node delete + edges to assigned + edges among remaining
        double h = 0.0;
        for (int ii : rem1) {
            h += _parameters.costs.node_del;
            for (int p = 0; p < n1; ++p) {
                int q = map12[p];
                if (q >= 0 && _source.edge(ii, p)) h += _parameters.costs.edge_del;
            }
        }
        for (int x = 0; x < (int)rem1.size(); ++x)
            for (int y = x+1; y < (int)rem1.size(); ++y)
                if (_source.edge(rem1[x], rem1[y])) h += _parameters.costs.edge_del;
        return h;
    }

    int k = std::max(r1, r2);
    std::vector<std::vector<double>> M(k, std::vector<double>(k, 0.0));

    // Fill real rows/cols
    for (int ri = 0; ri < k; ++ri) {
        bool row_is_node = (ri < r1);
        int i = row_is_node ? rem1[ri] : -1; // -1 denotes ε-row
        for (int cj = 0; cj < k; ++cj) {
            bool col_is_node = (cj < r2);
            int j = col_is_node ? rem2[cj] : -1; // -1 denotes ε-col
            double c = 0.0;
            if (row_is_node && col_is_node) {
                // i -> j substitution + edges to already-assigned
                c += nodeSubCost(i, j);
                // edges to assigned nodes:
                for (int p = 0; p < _source.nodes(); ++p) {
                    int q = map12[p];
                    if (q >= 0) {
                        bool e1 = _source.edge(i, p);
                        bool e2 = _target.edge(j, q);
                        c += edgeMismatchCost(e1, e2);
                    }
                }
            } else if (row_is_node && !col_is_node) {
                // i -> ε (deletion)
                c += _parameters.costs.node_del;
                for (int p = 0; p < _source.nodes(); ++p) {
                    int q = map12[p];
                    if (q >= 0 && _source.edge(i, p)) c += _parameters.costs.edge_del;
                }
            } else if (!row_is_node && col_is_node) {
                // ε -> j (insertion)
                c += _parameters.costs.node_ins;
                for (int p = 0; p < _source.nodes(); ++p) {
                    int q = map12[p];
                    if (q >= 0 && _target.edge(j, q)) c += _parameters.costs.edge_ins;
                }
            } else {
                // ε -> ε
                c = 0.0;
            }
            M[ri][cj] = std::max(0.0, c); // ensure non-negative
        }
    }

    auto res = hungarian_min(M);
    return res.first;
}

// Heuristic wrapper for a state
inline double GEDApproximation::heuristic_for_state(const State& s) {
  return heuristic_LSAP(s.map12, s.used2, s.next_i);
}

// Incremental edge cost for assigning G1 node i -> G2 node j (or j==-1 for deletion)
// against already-assigned nodes only (so we don't double-count).
inline double GEDApproximation::incEdgeCost_assign_i_to_j(const std::vector<int>& map12, int i, int j) const {
  double cost = 0.0;
  const int n1 = _source.nodes();
  for (int p = 0; p < n1; ++p) {
    int q = map12[p];
    if (p == i || q == -2) continue; // shouldn't happen
    if (q >= 0) { // p already mapped to q
      bool e1 = _source.edge(i, p);
      bool e2 = (j >= 0) ? _target.edge(j, q) : false;
      cost += edgeMismatchCost(e1, e2);
    }
    else if (q == -1) { // p was deleted
      bool e1 = _source.edge(i, p);
      if (e1) cost += _parameters.costs.edge_del;
    }
  }
  return cost;
}
