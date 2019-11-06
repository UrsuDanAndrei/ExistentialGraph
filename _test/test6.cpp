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
        "(A, [A, B])",
        "([A], [[A], B])",
        "([A, B], [[A, B], C])",
        "(A, B, C, D, [A, [B, C], [D, [A, [B]]]])"
    };

    std::vector<AEGraph> inputs;
    for (auto s : input_strs) {
        inputs.emplace_back(AEGraph(s));
    }

    std::vector<std::vector<int>> deiterations = {
        {0, 0},
        {1, 0},
        {1, 0},
        {0, 1, 1}
    };

    std::vector<AEGraph> refs = {
        AEGraph("([B], A)"),
        AEGraph("([A], [B])"),
        AEGraph("([A, B], [C])"),
        AEGraph("([[B, C], [[[B], A]], A], A, B, C, D)")
    };

    assert(inputs.size() == refs.size());

    std::cerr << "==================== Test 6 ===================\n";
    std::cerr << "Testing deiterate()...\n";
    size_t len = inputs.size();
    unsigned int total = len;
    for (size_t i = 0; i < len; i++) {
        auto out = inputs[i].deiterate(deiterations[i]);
        out.sort();
        refs[i].sort();

        if (out != refs[i]) {
            total--;
            std::cerr << "Wrong answer for input number " << i+1 << std::endl;

            std::cerr << "Graph: " << input_strs[i] << std::endl;
            for (auto e : deiterations[i])
                std::cerr << e << " ";
            std::cerr << std::endl;
            std::cerr << "Expected: " << refs[i].repr() << std::endl;
            std::cerr << "Got: " << out.repr() << std::endl << std::endl;
        }
    }

    if (total == len) {
        std::cerr << "passed: " << total * 2.5 << "/" << len * 2.5 << std::endl;
        std::cout << total * 2.5 << std::endl;
    } else {
        std::cerr << "failed: " << total * 2.5 << "/" << len * 2.5 << std::endl;
        std::cout << 0 << std::endl;
    }
    return 0;
}
