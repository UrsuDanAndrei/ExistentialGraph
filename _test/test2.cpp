// Copyright 2019 Luca Istrate, Danut Matei
#include <iostream>
#include <utility>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string>
#include "../aegraph.h"

int main() {
    std::vector<std::string> input_strs {
        "([[A]], [[P], B])",
        "([[[A]]], B)",
        "([A, [B], [[C]]])",
        "(A, [[B]], [C])",
        "(S, [[P]], [A, [B], [[C, D]]])"
    };

    std::vector<AEGraph> inputs;
    for (auto s : input_strs) {
        inputs.emplace_back(AEGraph(s));
    }

    std::vector<AEGraph> refs = {
        AEGraph("(A, [[P], B])"),
        AEGraph("([A], B)"),
        AEGraph("([A, [B], C])"),
        AEGraph("(A, B, [C])"),
        AEGraph("(S, [[P]], [A, [B], C, D])")
    };

    std::vector<std::vector<int>> cuts {
        {0},
        {0, 0},
        {0, 1},
        {1},
        {0, 1}
    };

    assert(inputs.size() == refs.size());

    std::cerr << "==================== Test 2 ===================\n";
    std::cerr << "Testing double_cut()...\n";
    size_t len = inputs.size();
    unsigned int total = len * 2;
    for (size_t i = 0; i < len; i++) {
        auto res = inputs[i].double_cut(cuts[i]);
        res.sort();
        refs[i].sort();
        if (res != refs[i]) {
            std::cerr << "Wrong answer for input number " << i+1 << std::endl;
            std::cerr << "Graph: " << input_strs[i] << std::endl;
            std::cerr << "Cut: ";
            for (auto e : cuts[i])
                std::cerr << e << " ";
            std::cerr << std::endl;

            std::cerr << "Expected: " << std::endl;
            std::cerr << refs[i].repr() << std::endl;
            std::cerr << "Got: " << std::endl;
            std::cerr << res.repr() << std::endl << std::endl;

            total -= 2;
        }
    }

    if (total == len * 2) {
        std::cerr << "passed: " << total << "/" << len * 2<< std::endl;
        std::cout << total << std::endl;
    } else {
        std::cerr << "failed: " << total << "/" << len * 2 << std::endl;
        std::cout << 0 << std::endl;
    }
    return 0;
}
