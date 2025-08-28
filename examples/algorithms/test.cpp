//
// Created by florian on 28.08.25.
//

#include <bits/stdc++.h>
using namespace std;

// ---------------- Graph structure ----------------
struct Graph {
    int n;
    vector<int> label;
    vector<vector<int>> adj;
    Graph(int n=0): n(n), label(n), adj(n, vector<int>(n,0)) {}
    void add_edge(int u,int v){ adj[u][v]=1; adj[v][u]=1; }
};

// ---------------- Cost model ----------------
struct CostModel {
    double node_ins=1, node_del=1, node_sub=1;
    double edge_ins=1, edge_del=1;
    double nodeSub(int a,int b) const {
        return (a==b)?0:node_sub;
    }
};

// ---------------- Node of search tree ----------------
struct Node {
    vector<pair<int,int>> mapping; // (u->v) or (u->-1 for deletion), (-1->v) insertion
    double g=0, lb=0;
    Node* parent=nullptr;
};

struct DF_GED {
    Graph g1,g2;
    CostModel C;
    vector<vector<double>> Cv, Ce;
    vector<int> sortedV1;
    double UB=1e18;
    vector<pair<int,int>> bestPath;

    DF_GED(const Graph& a,const Graph& b,const CostModel& cm): g1(a), g2(b), C(cm) {
        buildCostMatrices();
        sortVertices();
    }

    void buildCostMatrices(){
        int n=g1.n, m=g2.n;
        Cv.assign(n+2, vector<double>(m+2,1e9));
        // fill Cv with substitution, insertion, deletion costs
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++) Cv[i][j]=C.nodeSub(g1.label[i],g2.label[j]);
            Cv[i][m]=C.node_del; // u->ε
            Cv[n][j]=C.node_ins; // ε->v
        }
        // Ce analogous (here left simplified, can be refined)
        Ce.assign(n+2, vector<double>(m+2,0));
    }

    void sortVertices(){
        // Simple: sort by degree descending
        sortedV1.resize(g1.n);
        iota(sortedV1.begin(), sortedV1.end(),0);
        sort(sortedV1.begin(), sortedV1.end(),[&](int a,int b){
            int da=accumulate(g1.adj[a].begin(),g1.adj[a].end(),0);
            int db=accumulate(g1.adj[b].begin(),g1.adj[b].end(),0);
            return da>db;
        });
    }

    // Lower bound using assignment (placeholder, could plug Hungarian)
    double lowerBound(const Node& p){
        // For simplicity, return 0 (admissible but weak).
        return 0;
    }

    void dfs(Node* cur, vector<int>& used2, int depth){
        if(depth==g1.n){
            // Insert all remaining unmatched v ∈ V2
            Node complete=*cur;
            for(int v=0; v<g2.n; v++) if(!used2[v]) {
                complete.mapping.push_back({-1,v});
                complete.g+=C.node_ins;
            }
            if(complete.g<UB){ UB=complete.g; bestPath=complete.mapping; }
            return;
        }
        int u=sortedV1[depth];
        // try match u with v
        for(int v=0; v<g2.n; v++) if(!used2[v]){
            Node child=*cur;
            child.mapping.push_back({u,v});
            child.g+=Cv[u][v];
            child.lb=lowerBound(child);
            if(child.g+child.lb<UB){
                used2[v]=1;
                dfs(&child,used2,depth+1);
                used2[v]=0;
            }
        }
        // try delete u
        Node del=*cur;
        del.mapping.push_back({u,-1});
        del.g+=C.node_del;
        del.lb=lowerBound(del);
        if(del.g+del.lb<UB) dfs(&del,used2,depth+1);
    }

    double run(){
        UB=1e18;
        Node root;
        vector<int> used2(g2.n,0);
        dfs(&root,used2,0);
        return UB;
    }
};

// ---------------- Example ----------------
int main(){
    Graph g1(2), g2(2);
    g1.label={0,1}; g2.label={0,2};
    g1.add_edge(0,1); g2.add_edge(0,1);
    CostModel C;
    DF_GED solver(g1,g2,C);
    double dist=solver.run();
    cout<<"GED = "<<dist<<"\n";
    return 0;
}
