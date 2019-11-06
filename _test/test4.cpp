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
        "(A, B, C, D, [[E]])",
        "([P])",
        "(A, B)",
        "([[A, B]])",
        "([A], [B])",
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

    std::vector<std::vector<int>> erasures = {
        {3},
        {0},
        {0},
        {0, 0, 0},
        {0},
        {1},
        {1},
        {0},
        {2},
        {1}
    };

    std::vector<AEGraph> refs = {
        AEGraph("([[E]], A, B, D)"),
        AEGraph("()"),
        AEGraph("(B)"),
        AEGraph("([[B]])"),
        AEGraph("([B])"),
        AEGraph("([[A]])"),
        AEGraph("([[[A]]])"),
        AEGraph("()"),
        AEGraph("([C], [[B]])"),
        AEGraph("([[B], [[C, D]], A], S)")
    };

    assert(inputs.size() == refs.size());

    std::cerr << "==================== Test 4 ===================\n";
    std::cerr << "Testing erase()...\n";
    size_t len = inputs.size();
    unsigned int total = len;
    for (size_t i = 0; i < len; i++) {
        auto out = inputs[i].erase(erasures[i]);
        out.sort();
        refs[i].sort();

        if (out != refs[i]) {
            total--;
            std::cerr << "Wrong answer for input number " << i+1 << std::endl;

            std::cerr << "Graph: " << input_strs[i] << std::endl;
            for (auto e : erasures[i])
                std::cerr << e << " ";
            std::cerr << std::endl;
            std::cerr << "Expected: " << refs[i].repr() << std::endl;
            std::cerr << "Got: " << out.repr() << std::endl << std::endl;
        }
    }

    if (total == len) {
        std::cerr << "passed: " << total << "/" << len << std::endl;
        std::cout << total << std::endl;
    } else {
        std::cerr << "failed: " << total << "/" << len << std::endl;
        std::cout << 0 << std::endl;
    }
    return 0;
}
