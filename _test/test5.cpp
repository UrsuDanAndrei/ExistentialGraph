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
        "([A], [A, B])",
        "([A, B], [[A, B], C])",
        "(A, B, C, D, [A, [B, C], [D, [A, [B]]]])"
    };

    std::vector<AEGraph> inputs;
    for (auto s : input_strs) {
        inputs.emplace_back(AEGraph(s));
    }

    std::vector<std::vector<std::vector<int>>> refs = {
        {{0, 0}},
        {},
        {{1, 0}},
        {{0, 0, 0}, {0, 0, 1}, {0, 1, 0, 1}, {0, 1, 1}, {0, 2}}
    };

    assert(inputs.size() == refs.size());

    std::cerr << "==================== Test 5 ===================\n";
    std::cerr << "Testing possible_deiterations()...\n";
    size_t len = inputs.size();
    unsigned int total = len * 5;
    for (size_t i = 0; i < len; i++) {
        auto out = inputs[i].possible_deiterations();

        std::sort(out.begin(), out.end());
        std::sort(refs[i].begin(), refs[i].end());

        if (out != refs[i]) {
            total -= 5;
            std::cerr << "Wrong answer for input number " << i+1 << std::endl;
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

    if (total == len * 5) {
        std::cerr << "passed: " << total << "/" << len * 5 << std::endl;
        std::cout << total << std::endl;
    } else {
        std::cerr << "failed: " << total << "/" << len * 5 << std::endl;
        std::cout << 0 << std::endl;
    }
    return 0;
}
