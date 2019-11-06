// Copyright 2019 Luca Istrate, Danut Matei,
// Dan-Andrei Ursu 314CA, Adrian-Tudor Dumitrescu 314CA
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include "./aegraph.h"

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
	// returns all the paths in the tree that lead to possible double cuts
    std::vector<std::vector<int>> pdc;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
    	if (subgraphs[i].num_subgraphs() == 1 && subgraphs[i].size() == 1) {
    		pdc.push_back({i});
    	}

    	auto r = subgraphs[i].possible_double_cuts();
    	for (auto& v : r) {
    		v.insert(v.begin(), i);
    	}

    	// adding to the pdc of this graph the pdc of his subgraph i
    	copy(r.begin(), r.end(), back_inserter(pdc));
    }

    return pdc;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
	// returns an AEGraph with the double cut placed at the path <where> removed
    AEGraph g = *this;
    AEGraph *pg = &g;

    // reaching the end of the path
	int len_where = where.size();
	for (int i = 0; i < len_where - 1; i++) {
		pg = &pg->subgraphs[where[i]];
	}

	// preserving the subgraphs and atoms of the double cut
	int index = where[len_where - 1];
	copy(pg->subgraphs[index].subgraphs[0].subgraphs.begin(),
		pg->subgraphs[index].subgraphs[0].subgraphs.end(),
		back_inserter(pg->subgraphs));

	copy(pg->subgraphs[index].subgraphs[0].atoms.begin(),
		pg->subgraphs[index].subgraphs[0].atoms.end(),
		back_inserter(pg->atoms));

	// removing the double cut
	pg->subgraphs.erase(pg->subgraphs.begin() + index);

    return g;
}


std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
	// returns all the paths in the tree that lead to possible erasures
	std::vector<std::vector<int>> pe;

    int len_subgraphs = subgraphs.size();
    for (int i = 0; i < len_subgraphs; i++) {
		auto r = subgraphs[i].possible_erasures(level + 1);
		for (auto& v : r) {
			v.insert(v.begin(), i);
		}

		// adding to the pe of this graph the pe of his subgraph i
		copy(r.begin(), r.end(), back_inserter(pe));
	}

	// adding to the pe all the atoms and subgraphs on even levels
	if ((level + 2) % 2 == 1) {
		if (this->size() == 1 && !is_SA) {
			return pe;
		}

		for (int i = 0; i < len_subgraphs; i++) {
			pe.push_back({i});
		}

		int len_atoms = atoms.size();
		for (int i  = 0; i < len_atoms; i++) {
			pe.push_back({len_subgraphs + i});
		}
	}

	return pe;
}


AEGraph AEGraph::erase(std::vector<int> where) const {
    // returns an AEGraph with the erasures placed at the path <where>
    AEGraph g = *this;
    AEGraph *pg = &g;

    // reaching the end of the path
	int len_where = where.size();
	for (int i = 0; i < len_where - 1; i++) {
		pg = &pg->subgraphs[where[i]];
	}

	int index = where[len_where - 1];
	// checking if index is atom or subgraph to erase
	if (index >= pg->num_subgraphs()) {
		pg->atoms.erase(pg->atoms.begin() + index - pg->num_subgraphs());
	} else {
		pg->subgraphs.erase(pg->subgraphs.begin() + index);
	}

	return g;
}


std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
	// returns all the paths in the tree that lead to possible deiterations
    std::vector<std::vector<int>> pd;

    int len_subgraphs = subgraphs.size();
    for (int i = 0; i < len_subgraphs; i++) {
    	auto r = get_paths_to(subgraphs[i]);
    	auto index = r.begin();
    	int nr = 0;

    	for (auto it = r.begin(); it != r.end(); it++) {
    		if (it->size() == 1) {
    			index = it;
    			nr++;
    		}
    	}

    	// if subgraph i is the only subgraph alike it is not a pd
    	if (nr == 1) {
    		r.erase(index);
    	}

    	// adding as pd all subgraphs that are alike subghraph i
    	copy(r.begin(), r.end(), back_inserter(pd));

    	r = subgraphs[i].possible_deiterations();

    	for (auto& v : r) {
    		v.insert(v.begin(), i);
    	}

    	// adding to the pd of this graph the pd of his subgraph i
    	copy(r.begin(), r.end(), back_inserter(pd));
    }

    int len_atoms = atoms.size();
    for (int i = 0; i < len_atoms; i++) {
    	auto r = get_paths_to(atoms[i]);
    	auto index = r.begin();
    	int nr = 0;

    	for (auto it = r.begin(); it != r.end(); it++) {
    		if (it->size() == 1) {
    			index = it;
    			nr++;
    		}
    	}

    	// if atom i is the only atom alike it is not a pd
    	if (nr == 1) {
    		r.erase(index);
    	}

    	// adding as pd all atoms that are alike atom i
    	copy(r.begin(), r.end(), back_inserter(pd));
    }

    // leaving only unique pds
    for (auto it = pd.begin(); it != pd.end(); it++) {
    	pd.erase(std::remove(it + 1, pd.end(), *it), pd.end());
    }

    return pd;
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    // returns an AEGraph with the deiteration placed at the path <where>
    AEGraph g = *this;
    AEGraph *pg = &g;

    // reaching the end of the path
	int len_where = where.size();
	for (int i = 0; i < len_where - 1; i++) {
		pg = &pg->subgraphs[where[i]];
	}

	int index = where[len_where - 1];
	// checking if index is atom or subgraph to deiterate
	if (index >= pg->num_subgraphs()) {
		pg->atoms.erase(pg->atoms.begin() + index - pg->num_subgraphs());
	} else {
		pg->subgraphs.erase(pg->subgraphs.begin() + index);
	}

	return g;
}

