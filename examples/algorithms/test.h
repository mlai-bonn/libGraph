// DF_GED_Algorithm2_with_ops.cpp
// Depth-First GED (DF-GED, Algorithm 2) with explicit operation path.
// C++20 version with correct feasible UB and complete edge accounting.

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// ----------------------------- Graph -----------------------------
struct Graph {
    int n;
    bool directed{false};
    vector<int> label;                // node labels (use -1 if unlabeled)
    vector<vector<uint8_t>> adj;      // 0/1 adjacency

    Graph(int n_=0, bool directed_=false): n(n_), directed(directed_), label(n_, -1) {
        adj.assign(n, vector<uint8_t>(n, 0));
    }
    void add_edge(int u, int v) {
        adj[u][v] = 1;
        if (!directed) adj[v][u] = 1;
    }
    int degree(int u) const {
        int d = 0; for (int v=0; v<n; ++v) d += adj[u][v]; return d;
    }
};

// ----------------------------- Cost model -----------------------------
struct Cost {
    double c_v_ins = 1.0, c_v_del = 1.0;
    double c_e_ins = 1.0, c_e_del = 1.0;
    // Node substitution (override if you have attributes)
    double c_v_sub(int la, int lb) const { return (la == lb) ? 0.0 : 1.0; }
};

// ----------------------------- Hungarian (LSAP) -----------------------------
struct LSAPResult {
    vector<int> col_of_row; // assigned column j for each row i
    vector<int> row_of_col; // assigned row i for each column j
    double cost;
};

// O(n^3) min-cost assignment (Hungarian). Returns full assignment.
static LSAPResult MunkresSolve(const vector<vector<double>>& a) {
    const double INF = 1e100;
    int n = static_cast<int>(a.size());
    vector<double> u(n+1), v(n+1);
    vector<int> p(n+1), way(n+1);
    for (int i=1; i<=n; ++i) {
        p[0] = i;
        int j0 = 0;
        vector<double> minv(n+1, INF);
        vector<char> used(n+1, false);
        do {
            used[j0] = true;
            int i0 = p[j0], j1 = 0;
            double delta = INF;
            for (int j=1; j<=n; ++j) if (!used[j]) {
                double cur = a[i0-1][j-1] - u[i0] - v[j];
                if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (int j=0; j<=n; ++j) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        do { int j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
    }
    double value = 0.0;
    vector<int> col_of_row(n, -1), row_of_col(n, -1);
    for (int j=1; j<=n; ++j) {
        int i = p[j]-1, jj = j-1;
        value += a[i][jj];
        col_of_row[i] = jj;
        row_of_col[jj] = i;
    }
    return {col_of_row, row_of_col, value};
}

// ----------------------------- Operation recording -----------------------------
enum class OpType { SubstituteV, DeleteV, InsertV, DeleteE, InsertE };

struct Operation {
    OpType type;
    int a, b; // meaning depends on type:
              // SubstituteV: a=u in G1, b=v in G2
              // DeleteV:    a=u in G1, b=-1
              // InsertV:    a=-1, b=v in G2
              // DeleteE:    a=u, b=x endpoints in G1
              // InsertE:    a=v, b=y endpoints in G2
};

static string op_to_string(const Operation& op) {
    switch (op.type) {
        case OpType::SubstituteV: return "substitute " + to_string(op.a) + "->" + to_string(op.b);
        case OpType::DeleteV:     return "deleteV "   + to_string(op.a);
        case OpType::InsertV:     return "insertV "   + to_string(op.b);
        case OpType::DeleteE:     return "deleteE ("  + to_string(op.a) + "," + to_string(op.b) + ") (G1)";
        case OpType::InsertE:     return "insertE ("  + to_string(op.a) + "," + to_string(op.b) + ") (G2)";
    }
    return {};
}

// ----------------------------- Edit Path (node in search tree) -----------------------------
struct EditPath {
    // mapping[u] = v >=0: u->v ; = -1: u->ε ; = -2: undecided
    vector<int> mapping;
    // inv2[v] = u >=0: v is image of u ; = -1: free (unmapped) ; = -2: inserted
    vector<int> inv2;

    int depth{0};        // number of decided vertices of V1 (by sorted_V1 order)
    double g{0.0}, lb{0.0};
    int parent_id{-1};
    vector<Operation> ops; // operation sequence to reach this node
};

// ----------------------------- DF-GED solver (Algorithm 2) -----------------------------
struct DF_GED {
    // Inputs
    Graph g1, g2; Cost C;

    // Preprocessing artifacts
    vector<vector<double>> Cv, Ce;   // (Cv used for LB; Ce placeholder)
    vector<int> sorted_V1;

    // Global search state
    double UB{numeric_limits<double>::infinity()};
    vector<Operation> Best_Edit_Path; // explicit ops
    vector<int> OPEN;                 // indices of candidate nodes
    vector<EditPath> pool;            // pool of nodes

    DF_GED(const Graph& G1, const Graph& G2, const Cost& cost): g1(G1), g2(G2), C(cost) {
        generate_Cv_Ce_and_sorted_V1();
    }

    // Preprocessing: generate Cv, Ce, and sorted_V1
    void generate_Cv_Ce_and_sorted_V1() {
        const double INF = numeric_limits<double>::infinity();
        int n = g1.n, m = g2.n;

        Cv.assign(n+2, vector<double>(m+2, INF));
        for (int i=0;i<n;++i)
            for (int j=0;j<m;++j)
                Cv[i][j] = C.c_v_sub(g1.label[i], g2.label[j]);
        for (int i=0;i<n;++i) Cv[i][m] = C.c_v_del; // u->ε
        for (int j=0;j<m;++j) Cv[n][j] = C.c_v_ins; // ε->v
        Ce.assign(n+2, vector<double>(m+2, INF));   // placeholder (edge LB not used here)

        sorted_V1.resize(n);
        iota(sorted_V1.begin(), sorted_V1.end(), 0);
        sort(sorted_V1.begin(), sorted_V1.end(), [&](int a, int b){
            int da = g1.degree(a), db = g1.degree(b);
            if (da != db) return da > db;
            return a < b;
        });
    }

    // Pending sets
    vector<int> pending_vertices1(const EditPath& P) const {
        vector<int> out; out.reserve(g1.n);
        for (int u=0; u<g1.n; ++u) if (P.mapping[u] == -2) out.push_back(u);
        return out;
    }
    vector<int> pending_vertices2(const EditPath& P) const {
        vector<int> out; out.reserve(g2.n);
        for (int v=0; v<g2.n; ++v) if (P.inv2[v] == -1) out.push_back(v);
        return out;
    }

    // Lower bound lb(p): node-only LSAP over remaining vertices (admissible)
    double lb_of(const EditPath& P) const {
        vector<int> U1, U2;
        for (int u=0; u<g1.n; ++u) if (P.mapping[u] == -2) U1.push_back(u);
        for (int v=0; v<g2.n; ++v) if (P.inv2[v] == -1) U2.push_back(v);

        int n1 = static_cast<int>(U1.size()), n2 = static_cast<int>(U2.size());
        int L = max(n1, n2);
        if (L == 0) return 0.0;

        vector<vector<double>> A(L, vector<double>(L, 0.0));
        for (int i=0; i<L; ++i) for (int j=0; j<L; ++j) {
            if (i < n1 && j < n2) {
                int u = U1[i], v = U2[j];
                A[i][j] = C.c_v_sub(g1.label[u], g2.label[v]);
            } else if (i < n1 && j >= n2) {
                A[i][j] = C.c_v_del; // extra rows -> deletions
            } else if (i >= n1 && j < n2) {
                A[i][j] = C.c_v_ins; // extra cols -> insertions
            } else {
                A[i][j] = 0.0;
            }
        }
        return MunkresSolve(A).cost;
    }

    // Incremental cost + operations when deciding u -> v (v >= 0) or u -> ε (v == -1).
    double g_of_extend_and_ops(const EditPath& P, int u, int v, vector<Operation>& newOps) const {
        double dc = 0.0;
        if (v == -1) {
            // Delete vertex u. Charge edges to already-decided mapped OR deleted neighbors.
            dc += C.c_v_del;
            newOps.push_back({OpType::DeleteV, u, -1});
            for (int x=0; x<g1.n; ++x) {
                if (!g1.adj[u][x]) continue;
                int y = P.mapping[x];
                if (y >= 0 || y == -1) { // neighbor already mapped or deleted
                    dc += C.c_e_del;
                    newOps.push_back({OpType::DeleteE, u, x});
                }
                // if y == -2 (undecided), we charge later when x gets decided
            }
            return P.g + dc;
        }

        // Substitute u -> v
        dc += C.c_v_sub(g1.label[u], g2.label[v]);
        newOps.push_back({OpType::SubstituteV, u, v});

        for (int x=0; x<g1.n; ++x) if (x != u) {
            int y = P.mapping[x];
            if (y == -2) continue; // undecided counterpart: handle later
            if (y == -1) {
                // counterpart x was deleted earlier; any edge (u,x) in G1 must be deleted now
                if (g1.adj[u][x]) {
                    dc += C.c_e_del;
                    newOps.push_back({OpType::DeleteE, u, x});
                }
            } else { // y >= 0 (mapped)
                uint8_t e1 = g1.adj[u][x];
                uint8_t e2 = g2.adj[v][y];
                if (e1 && !e2) {
                    dc += C.c_e_del;
                    newOps.push_back({OpType::DeleteE, u, x});
                } else if (!e1 && e2) {
                    dc += C.c_e_ins;
                    newOps.push_back({OpType::InsertE, v, y});
                }
            }
        }
        return P.g + dc;
    }

    // bestChild(parenttmp): child in OPEN with minimal g+lb whose parent == parenttmp
    int bestChild(int parent_id) const {
        double bestF = numeric_limits<double>::infinity();
        int bestId = -1;
        for (int id : OPEN) {
            const EditPath& Q = pool[id];
            if (Q.parent_id == parent_id) {
                double f = Q.g + Q.lb;
                if (f < bestF) { bestF = f; bestId = id; }
            }
        }
        return bestId; // -1 denotes φ
    }

    // backtrack(parenttmp): parent of a node
    int backtrack(int parent_id) const {
        if (parent_id < 0) return -1;
        return pool[parent_id].parent_id;
    }

    void remove_from_OPEN(int id) {
        auto it = find(OPEN.begin(), OPEN.end(), id);
        if (it != OPEN.end()) OPEN.erase(it);
    }

    void maybe_push_to_OPEN(EditPath&& child) {
        if (child.g + child.lb < UB) {
            int id = static_cast<int>(pool.size());
            pool.push_back(std::move(child));
            OPEN.push_back(id);
        }
    }

    // Expand a node like Algorithm 2 (lines 18–26).
    void expand(int pmin_id) {
        EditPath P = pool[pmin_id]; // snapshot for branching

        auto pend1 = pending_vertices1(P);
        if (!pend1.empty()) {
            int u_next = sorted_V1[P.depth];

            // Branch over uk+1 -> w for each free w in V2
            auto pend2 = pending_vertices2(P);
            for (int w : pend2) {
                EditPath Q = P;
                Q.mapping[u_next] = w;
                Q.inv2[w] = u_next;
                Q.depth = P.depth + 1;
                vector<Operation> deltaOps;
                Q.g = g_of_extend_and_ops(P, u_next, w, deltaOps);
                Q.ops.insert(Q.ops.end(), deltaOps.begin(), deltaOps.end());
                Q.lb = lb_of(Q);
                Q.parent_id = pmin_id;
                maybe_push_to_OPEN(std::move(Q));
            }
            // Branch uk+1 -> ε
            {
                EditPath Q = P;
                Q.mapping[u_next] = -1;
                Q.depth = P.depth + 1;
                vector<Operation> deltaOps;
                Q.g = g_of_extend_and_ops(P, u_next, -1, deltaOps);
                Q.ops.insert(Q.ops.end(), deltaOps.begin(), deltaOps.end());
                Q.lb = lb_of(Q);
                Q.parent_id = pmin_id;
                maybe_push_to_OPEN(std::move(Q));
            }
        } else {
            // Completion: insert remaining vertices ε -> w and their incident edges
            EditPath Q = P;
            auto pend2 = pending_vertices2(Q); // vertices to insert in V2
            // Insert remaining nodes
            for (int w : pend2) {
                Q.ops.push_back({OpType::InsertV, -1, w});
                Q.g += C.c_v_ins;
                Q.inv2[w] = -2; // mark as inserted
            }
            // Insert edges between inserted nodes and mapped nodes
            // (a) inserted–mapped
            for (int w : pend2) {
                for (int y=0; y<g2.n; ++y) {
                    if (Q.inv2[y] >= 0 /*mapped*/ && g2.adj[w][y]) {
                        Q.ops.push_back({OpType::InsertE, w, y});
                        Q.g += C.c_e_ins;
                    }
                }
            }
            // (b) inserted–inserted
            for (int i=0; i<(int)pend2.size(); ++i)
                for (int j=i+1; j<(int)pend2.size(); ++j)
                    if (g2.adj[pend2[i]][pend2[j]]) {
                        Q.ops.push_back({OpType::InsertE, pend2[i], pend2[j]});
                        Q.g += C.c_e_ins;
                    }

            // Update UB & Best_Edit_Path if improved
            if (Q.g < UB) {
                UB = Q.g;
                Best_Edit_Path = Q.ops; // full explicit operations
            }
        }
    }

    // Seed OPEN with children of the root for u1 in sorted_V1 (lines 4–5)
    void seed_OPEN_with_root_children() {
        EditPath root;
        root.mapping.assign(g1.n, -2);
        root.inv2.assign(g2.n, -1);
        root.depth = 0;
        root.g = 0.0;
        root.lb = lb_of(root);
        root.parent_id = -1;
        int root_id = static_cast<int>(pool.size());
        pool.push_back(root);

        if (g1.n == 0) return; // degenerate

        int u1 = sorted_V1[0];
        // u1 -> w
        for (int w=0; w<g2.n; ++w) {
            EditPath P = root;
            P.mapping[u1] = w;
            P.inv2[w] = u1;
            P.depth = 1;
            vector<Operation> deltaOps;
            P.g = g_of_extend_and_ops(root, u1, w, deltaOps);
            P.ops.insert(P.ops.end(), deltaOps.begin(), deltaOps.end());
            P.lb = lb_of(P);
            P.parent_id = root_id;
            maybe_push_to_OPEN(std::move(P));
        }
        // u1 -> ε
        {
            EditPath P = root;
            P.mapping[u1] = -1;
            P.depth = 1;
            vector<Operation> deltaOps;
            P.g = g_of_extend_and_ops(root, u1, -1, deltaOps);
            P.ops.insert(P.ops.end(), deltaOps.begin(), deltaOps.end());
            P.lb = lb_of(P);
            P.parent_id = root_id;
            maybe_push_to_OPEN(std::move(P));
        }
    }

    // Build a FEASIBLE initial UB from LSAP (node + all implied edge ops)
    pair<double, vector<Operation>>
    initial_upper_bound_from_LSAP() const {
        const double BIG = 1e9;
        int n = g1.n, m = g2.n, L = n + m;
        vector<vector<double>> A(L, vector<double>(L, BIG));

        // substitutions block
        for (int i=0;i<n;++i)
            for (int j=0;j<m;++j)
                A[i][j] = C.c_v_sub(g1.label[i], g2.label[j]);

        // deletions (rows i<n to cols m+i)
        for (int i=0;i<n;++i) A[i][m+i] = C.c_v_del;

        // insertions (rows n+j to cols j<m)
        for (int j=0;j<m;++j) A[n+j][j] = C.c_v_ins;

        auto res = MunkresSolve(A);
        const auto& col_of_row = res.col_of_row;
        const auto& row_of_col = res.row_of_col;

        vector<int> map1(n, -2); // u -> v or -1
        vector<int> inv2(m, -1); // v -> u or -1
        vector<int> inserted(m, 0), deleted(n, 0);
        vector<Operation> ops;
        double cost = 0.0;

        // decode assignment → vertex ops
        for (int i=0;i<n;++i) {
            int c = col_of_row[i];
            if (c < m) { // substitute u_i -> v_c
                map1[i] = c; inv2[c] = i;
                cost += C.c_v_sub(g1.label[i], g2.label[c]);
                ops.push_back({OpType::SubstituteV, i, c});
            } else {     // delete u_i
                map1[i] = -1; deleted[i] = 1;
                cost += C.c_v_del;
                ops.push_back({OpType::DeleteV, i, -1});
            }
        }
        for (int j=0;j<m;++j) {
            int r = row_of_col[j];
            if (r >= n) { // insert v_j
                inserted[j] = 1;
                cost += C.c_v_ins;
                ops.push_back({OpType::InsertV, -1, j});
            }
        }

        // edge ops:
        // (1) mapped-mapped pairs
        for (int u=0; u<n; ++u) if (map1[u] >= 0)
        for (int x=u+1; x<n; ++x) if (map1[x] >= 0) {
            int v = map1[u], y = map1[x];
            uint8_t e1 = g1.adj[u][x], e2 = g2.adj[v][y];
            if (e1 && !e2) { cost += C.c_e_del; ops.push_back({OpType::DeleteE, u, x}); }
            else if (!e1 && e2) { cost += C.c_e_ins; ops.push_back({OpType::InsertE, v, y}); }
        }
        // (2) deleted–mapped edges from G1
        for (int u=0; u<n; ++u) if (deleted[u])
        for (int x=0; x<n; ++x) if (map1[x] >= 0) {
            if (g1.adj[u][x]) { cost += C.c_e_del; ops.push_back({OpType::DeleteE, u, x}); }
        }
        // (3) inserted–mapped edges from G2
        for (int v=0; v<m; ++v) if (inserted[v])
        for (int y=0; y<m; ++y) if (inv2[y] >= 0) {
            if (g2.adj[v][y]) { cost += C.c_e_ins; ops.push_back({OpType::InsertE, v, y}); }
        }
        // (4) inserted–inserted edges
        for (int v=0; v<m; ++v) if (inserted[v])
        for (int y=v+1; y<m; ++y) if (inserted[y]) {
            if (g2.adj[v][y]) { cost += C.c_e_ins; ops.push_back({OpType::InsertE, v, y}); }
        }

        return {cost, ops};
    }

    // -------------------- Algorithm 2 main loop --------------------
    pair<double, vector<Operation>> run() {
        OPEN.clear();
        Best_Edit_Path.clear();
        pool.clear();

        // (line 3) UB from a FEASIBLE mapping (node+edge ops)
        {
            auto [ub0, ops0] = initial_upper_bound_from_LSAP();
            UB = ub0;
            Best_Edit_Path = std::move(ops0);
        }

        // (lines 4–5) seed OPEN with children of u1 ∈ sorted_V1
        seed_OPEN_with_root_children();

        int r = -1;            // parent(u1) for the root is -1
        int parenttmp = r;     // (line 6)

        while (true) {
            int pmin = bestChild(parenttmp);      // (line 8)

            // (lines 9–12)
            while (pmin == -1 && parenttmp != r) {
                parenttmp = backtrack(parenttmp);
                pmin = bestChild(parenttmp);
            }

            // (lines 13–15)
            if (pmin == -1 && parenttmp == r) {
                return {UB, Best_Edit_Path};
            }

            remove_from_OPEN(pmin);              // (line 16)

            EditPath& Pmin = pool[pmin];
            if (Pmin.g + Pmin.lb < UB) {         // (line 17)
                if (!pending_vertices1(Pmin).empty()) {
                    expand(pmin);                 // (lines 18–24)
                } else {
                    expand(pmin);                 // (lines 25–29)
                }
            }
            parenttmp = pmin;                     // (line 32)
        }
    }
};

// ----------------------------- Demo -----------------------------
static void test() {
    // Example: G1 triangle vs G2 path (all labels 0)
    Graph g1(3,false); g1.add_edge(0,1); g1.add_edge(1,2); g1.add_edge(0,2);
    Graph g2(3,false); g2.add_edge(0,1); g2.add_edge(1,2);
    g1.label = {0,0,0};
    g2.label = {0,0,0};

    Cost C; // unit costs

    DF_GED solver(g1, g2, C);
    auto [dist, ops] = solver.run();

    cout << fixed << setprecision(6);
    cout << "GED (UB) = " << dist << "\n";
    cout << "Edit Path (operations):\n";
    for (const auto& op : ops) {
        cout << "  - " << op_to_string(op) << "\n";
    }
}
