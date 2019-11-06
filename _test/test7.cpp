// Copyright 2019 Luca Istrate, Danut Matei
#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>
#include <string>
#include "../aegraph.h"

std::vector<std::pair<std::string, std::vector<int>>> reduce(AEGraph graph);

std::vector<std::pair<std::string, std::vector<int>>>
bactracking_step(AEGraph graph, std::string op) {
    if (op == "DE") {
        auto steps = graph.possible_deiterations();
        sort(steps.begin(), steps.end());
        for (auto &step : steps) {
            auto g = graph.deiterate(step);
            g.sort();
            auto r = reduce(g);
            if (!r.empty()) {
                r.insert(r.begin(), make_pair("DE", step));
                return r;
            }
        }
    } else if (op == "DC") {
        auto steps = graph.possible_double_cuts();
        sort(steps.begin(), steps.end());
        for (auto &step : steps) {
            auto g = graph.double_cut(step);
            g.sort();
            auto r = reduce(g);
            if (!r.empty()) {
                r.insert(r.begin(), make_pair("DC", step));
                return r;
            }
        }
    } else if (op == "E") {
        auto steps = graph.possible_erasures();
        sort(steps.begin(), steps.end());
        for (auto &step : steps) {
            auto g = graph.erase(step);
            g.sort();
            auto r = reduce(g);
            if (!r.empty()) {
                r.insert(r.begin(), make_pair("E", step));
                return r;
            }
        }
    }

    return {};
}

bool is_contradiction(AEGraph graph) {
    auto g = graph;

    if (g.num_atoms() != 1 || g.num_subgraphs() != 1)
        return false;

    if (g.subgraphs[0].num_subgraphs() != 0 || g.subgraphs[0].num_atoms() != 1)
        return false;

    return g.atoms[0] == g.subgraphs[0].atoms[0];
}

std::vector<std::pair<std::string, std::vector<int>>> reduce(AEGraph graph) {
    if (is_contradiction(graph))
        return {{"END", {}}};

    for (auto op : {"DC", "DE", "E"}) {
        auto r = bactracking_step(graph, op);
        if (!r.empty())
            return r;
    }

    return {};
}

AEGraph make_counterset(AEGraph premise, AEGraph conclusion) {
    std::string prem = premise.repr().substr(1, premise.repr().size() - 2);
    std::string conc = conclusion.repr().substr(1,
        conclusion.repr().size() - 2);
    return AEGraph("(" + prem + ", [" + conc + "])");
}

std::vector<std::pair<std::string, std::vector<int>>> steps_to(AEGraph premise,
    AEGraph conclusion) {
    auto graph = make_counterset(premise, conclusion);

    return reduce(graph);
}

int main() {
    std::cerr << "==================== Test 7 ===================\n";
    std::cerr << "Testing full proof...\n";
    std::vector<std::string> premises = {
        "(A)",
        "(A)",
        "(p, [p, [q]])",
        "(r, [p, [q]])",
        "(P, Q, R, [[A], [B]])"
    };

    std::vector<std::string> conclusions = {
        "([[A]])",
        "([A])",
        "(q)",
        "(q)",
        "(P, A, B)"
    };

    std::vector<int> has_solution = {1, 0, 1, 0, 0};

    size_t len = premises.size();
    size_t total = len;
    for (size_t i = 0; i < len; i++) {
        auto premise = AEGraph(premises[i]);
        auto conclusion = AEGraph(conclusions[i]);

        auto res = steps_to(premise, conclusion);
        if (has_solution[i]) {
            if (res.size() == 0) {
                total--;
                std::cerr << "Wrong answer for input number " << i+1
                    << std::endl;
                std::cerr << "Expected a solution" << std::endl;
            }
        } else {
            if (res.size() != 0) {
                total--;
                std::cerr << "Wrong answer for input number " << i+1
                    << std::endl;
                std::cerr << "Expected no solution" << std::endl;
            }
        }
    }

    if (total == len) {
        std::cerr << "passed: " << total * 2 << "/" << len * 2 << std::endl;
        std::cout << total * 2 << std::endl;
    } else {
        std::cerr << "failed: " << total * 2 << "/" << len * 2 << std::endl;
        std::cout << 0 << std::endl;
    }
    return 0;
}
