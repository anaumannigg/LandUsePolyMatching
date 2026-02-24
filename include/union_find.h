#ifndef TCPOLYGONMATCHING_UNION_FIND_H
#define TCPOLYGONMATCHING_UNION_FIND_H

#include <vector>

//minimalistic UnionFind datastructure
class UnionFind {
public:
    explicit UnionFind(size_t n)
        : parent(n), rank(n, 0)
    {
        for (size_t i = 0; i < n; ++i)
            parent[i] = i;
    }

    int find(int x) {
        if (parent[x] != x)
            parent[x] = find(parent[x]);  // path compression
        return parent[x];
    }

    void unite(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX == rootY)
            return;

        // union by rank
        if (rank[rootX] < rank[rootY]) {
            parent[rootX] = rootY;
        } else if (rank[rootX] > rank[rootY]) {
            parent[rootY] = rootX;
        } else {
            parent[rootY] = rootX;
            rank[rootX]++;
        }
    }

private:
    std::vector<int> parent;
    std::vector<int> rank;
};

#endif //TCPOLYGONMATCHING_UNION_FIND_H