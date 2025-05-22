//
// Created by junior on 2020/12/4.
//

#ifndef AVS_API_GRAPH_H
#define AVS_API_GRAPH_H

#include "dna.h"
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/static_graph.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <lemon/kruskal.h>

struct Weight {
    int off;
    int mismatch;

    Weight(int off, int mismatch) : off(off), mismatch(mismatch) {}

    explicit Weight() : off(0), mismatch(0) {}

    Weight(int value) : off(value), mismatch(value) {}

    friend bool operator<(const Weight &w1, const Weight &w2) {
        if (w1.off < w2.off) return true;
        else if (w1.off > w2.off) return false;
        else return w1.mismatch < w2.mismatch;
    }

    friend bool operator>(const Weight &w1, const Weight &w2) {
        if (w1.off > w2.off) return true;
        else if (w1.off < w2.off) return false;
        else return w1.mismatch > w2.mismatch;
    }

    friend bool operator==(const Weight &w1, const Weight &w2) {
        return w1.off == w2.off && w1.mismatch == w2.mismatch;
    }

    Weight &operator+=(const Weight &w) {
        off += w.off;
        mismatch += w.mismatch;
        return *this;
    }

    Weight &operator-=(const Weight &w) {
        off -= w.off;
        mismatch -= w.mismatch;
        return *this;
    }
};

struct Strand {
    u8 start_strand_id: 4, end_strand_id: 4; // e: start -> end

    explicit Strand() : start_strand_id(0), end_strand_id(0) {}

    Strand(u32 a, u32 b) : start_strand_id(a), end_strand_id(b) {}
};


struct Graph {
    using G = lemon::SmartDigraph;
    using WeightMap = G::ArcMap<Weight>;
    using StrandMap = G::ArcMap<Strand>;
    G g;
    WeightMap weight_map;
    StrandMap strand_map;

    explicit Graph(size_t read_count) : weight_map(g), strand_map(g) {
        g.reserveNode(read_count);
        for (size_t i = 0; i < read_count; i++) {
            g.addNode();
        }
    }
};

#endif //AVS_API_GRAPH_H
