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
        "(A)",
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

    std::vector<std::vector<std::vector<int>>> refs = {
        {},             // (A)
        {},             // ([P])
        {},             // (A, B)
        {{0}},          // ([[A, B]])
        {},             // ([A], [B])
        {{0}},          // ([[A]], [[P], B])
        {{0}, {0, 0}},  // ([[[A]]], B)
        {{0, 1}},       // ([A, [B], [[C]]])
        {{1}},          // (A, [[B]], [C])
        {{0, 1}, {1}}   // (S, [[P]], [A, [B], [[C, D]]])
    };

    assert(inputs.size() == refs.size());

    std::cerr << "==================== Test 1 ===================\n";
    std::cerr << "Testing possible_double_cuts()...\n";
    size_t len = inputs.size();
    unsigned int total = len;
    for (size_t i = 0; i < len; i++) {
        auto out = inputs[i].possible_double_cuts();
        sort(out.begin(), out.end());
        sort(refs[i].begin(), refs[i].end());
        if (out != refs[i]) {
            total--;
            std::cerr << "Wrong answer for input number " << i+1 << std::endl;
            std::cerr << "Graph: " << input_strs[i] << std::endl;

            std::cerr << "Expected:" << std::endl;
            for (auto v : refs[i]) {
                for (auto e : v) {
                    std::cerr << e << " ";
                }
                std::cerr << std::endl;
            }

            std::cerr << "Got:" << std::endl;
            for (auto v : out) {
                for (auto e : v) {
                    std::cerr << e << " ";
                }
                std::cerr << std::endl;
            }
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
