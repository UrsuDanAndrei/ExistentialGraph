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
        {{0}},                        // (A)
        {{0}},                        // ([P])
        {{0}, {1}},                   // (A, B)
        {{0}, {0, 0, 0}, {0, 0, 1}},  // ([[A, B]])
        {{0}, {1}},                   // ([A], [B])
        {{0}, {1}},                   // ([[A]], [[P], B])
        {{0}, {1}},                   // ([[[A]]], B)
        {{0}},                        // ([A, [B], [[C]]])
        {{0}, {1}, {2}},              // (A, [[B]], [C])
        {{0}, {1}, {2}}               // (S, [[P]], [A, [B], [[C, D]]])
    };

    assert(inputs.size() == refs.size());

    std::cerr << "==================== Test 3 ===================\n";
    std::cerr << "Testing possible_erasures()...\n";
    size_t len = inputs.size();
    unsigned int total = len;
    for (size_t i = 0; i < len; i++) {
        auto out = inputs[i].possible_erasures();
        std::sort(out.begin(), out.end());
        std::sort(refs[i].begin(), refs[i].end());

        if (out != refs[i]) {
            total--;
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

    if (total == len) {
        std::cerr << "passed: " << total << "/" << len << std::endl;
        std::cout << total << std::endl;
    } else {
        std::cerr << "failed: " << total << "/" << len << std::endl;
        std::cout << 0 << std::endl;
    }
    return 0;
}
