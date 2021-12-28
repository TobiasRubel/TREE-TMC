#include <chrono>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iomanip>
#include "heuristics/maxcut/burer2002.h"
#include "problem/instance.h"
#include "problem/max_cut_instance.h"

int verbose = 0, stats[16];

/*
    succinct names for unordered map and set of strings
    str<int>::map <=> std::unordered_map<std::string, int>
    str<void>::set <=> std::unordered_set<std::string>
*/

template <typename T> class str {
    public: 
        typedef std::unordered_map<std::string, T> map;
        typedef std::unordered_set<std::string> set;
};

/*
    basic n x n x 2 matrix processing functions
    a quartet graph is represented by two matrices G and B
    storing weights of good and bad edges, respectively
*/

class Matrix {
    public:
        template <typename T> static T ***new_mat(int size);
        template <typename T> static void delete_mat(T ***m, int size);
        template <typename T> static void display_mat(T ***m, int size, int k);
        template <typename T> static T diff_mat(T ***m1, T ***m2, int size, int k);
};

/*
    data structures for taxa and taxon relationship
    in the tree structure, each node represents a taxon
    inner nodes are artificial taxa while leaves are input taxa
    the node of X becomes the parent of the nodes of taxa in S
    when X represents the taxa in the set S
*/

class Taxa {
    public:
        class Node {
            friend class Taxa;
            public:
                Node(std::string label);
            private:
                Node *parent;
                std::string label;
                int index, size;
                bool is_single;
                double weight;
        };
        Taxa(str<void>::set labels, int weighting);
        Taxa(const Taxa &source);
        ~Taxa();
        std::string to_string();
        std::string info();
        int get_weighting();
        void update(str<void>::set labels, std::string pseudo);
        int leaf_size();
        int size();
        int single_size();
        int multiple_size();
        std::string index2label(int i);
        std::string index2leaflabel(int i);
        int label2index(std::string label);
        int label2key(std::string label);
        std::string key2label(int i);
        std::string get_label(std::string label);
        bool is_single(std::string label);
        double get_weight(std::string label);
    private: 
        std::vector<Node *> roots, leaves;
        str<Node *>::map label2node;
        int singles, weighting;
        Node *get_root(std::string label);
        Node *get_root(Node *node);
        double get_weight(Node *node);
};

/*
    data structures for gene trees and species trees
    rooted, binary, and multi-labeled
    each tree node is augmented with a vector and a matrix
    the vector stores the number of artificial taxa in the subtree
    the matrix stores the number of tuples forming quartets of certain topology
*/

class Tree {
    public: 
        class Node {
            friend class Tree;
            public:
                Node(const std::string &name);
                void new_states(int size);
                static std::string pseudonym();
                void delete_states();
                ~Node();
            private:
                static int pseudonyms;
                Node *left, *right, *parent;
                int size;
                double *leaves, **pairs, s1, s2;
                std::string label;
                static double get_pairs(double *leaves, double s1, double s2, int x, int y);
        };
        Tree(std::ifstream &fin, int execution, int weighting, int s0, int s1);
        Tree(const std::string &newick);
        std::string to_string();
        ~Tree();
        double ***build_graph(Taxa &subset);
        void append_labels(Node *root, str<void>::set &labels);
        void append_quartets(str<double>::map &quartets, Taxa &subset);
        static void append_quartet(str<double>::map &quartets, std::string quartet, double weight);
        static std::string join(std::string *labels);
        static std::string *split(const std::string &qt);
    private: 
        str<Node *>::map label2node;
        Node *root;
        void clear_states(Node *root);
        void build_states(Node *root, Taxa &subset);
        void build_depth(Node *root, int depth);
        double pairs(Node *subtree, int x, int y, bool comple);
        void single_pairs(Node *root, double s, int x);
        double multiple_pairs(Node *root, int x, int y);
        str<void>::set build_mat(Node *root, Taxa &subset, double ***mat);
        Node *build_tree(const std::string &newick);
        std::string display_tree(Node *root);
        Node *construct_stree(std::vector<Tree *> &input, Taxa &subset);
        Node *construct_stree_brute(str<double>::map &input, Taxa &subset);
        Node *construct_stree_check(str<double>::map &input, std::vector<Tree *> &input_, Taxa &subset);
        Node *reroot(Node *root, str<void>::set &visited);
        Node *reroot_stree(Node *root, const std::string &pseudo);
        Node *pseudo2node(Node *root, const std::string &pseudo);
        static std::string ordered(const std::string &a, const std::string &b);
};

/*
    data structures for quartet graphs
    the graph is built based on either quartets or gene trees
    rank-two ralaxation heuristics are used to compute max-cut 
*/

class Graph {
    public:
        Graph(std::vector<Tree*> &input, Taxa &subset);
        Graph(str<double>::map &input, Taxa &subset);
        double distance(Graph *g, int k);
        double get_cut(str<void>::set *A, str<void>::set *B);
        ~Graph();
    private:
        int size;
        str<int>::map label2index;
        std::vector<std::string> labels;
        double ***graph;
        double sdp_cut(double alpha, str<void>::set *A, str<void>::set *B);
};

template <typename T>
T ***Matrix::new_mat(int size) {
    T ***m = new T**[size];
    for (int i = 0; i < size; i ++) {
        m[i] = new T*[size];
        for (int j = 0; j < size; j ++) {
            m[i][j] = new T[2];
            m[i][j][0] = m[i][j][1] = 0;
        }
    }
    return m;
}

template <typename T>
void Matrix::delete_mat(T ***m, int size) {
    for (int i = 0; i < size; i ++) {
        for (int j = 0; j < size; j ++) {
            delete[] m[i][j];
        }
        delete[] m[i];
    }
    delete[] m;
}

template <typename T>
void Matrix::display_mat(T ***m, int size, int k) {
    if (k != 0 && k != 1) {
        for (int i = 0; i < size; i ++) {
            for (int j = 0; j < size; j ++) {
                std::cout << std::setw(8) << std::setprecision(6) << m[i][j][0] + m[i][j][1];
            }
            std::cout << std::endl;
        }
    }
    else {
        for (int i = 0; i < size; i ++) {
            for (int j = 0; j < size; j ++) {
                std::cout << std::setw(8) << std::setprecision(6) << m[i][j][k];
            }
            std::cout << std::endl;
        }
    }
}

template <typename T>
T Matrix::diff_mat(T ***m1, T ***m2, int size, int k) {
    T sum = 0;
    for (int i = 0; i < size; i ++) {
        for (int j = 0; j < size; j ++) {
            T delta = m1[i][j][k] - m2[i][j][k];
            if (delta < 0) delta = - delta;
            sum += delta;
        }
    }
    return sum;
}

Taxa::Node::Node(std::string label) {
    this->label = label;
    index = -1;
    parent = NULL;
}

Taxa::Taxa(str<void>::set labels, int weighting) {
    singles = 0;
    this->weighting = weighting;
    for (std::string label : labels) {
        Node *node = new Node(label);
        label2node[label] = node;
        roots.push_back(node);
        leaves.push_back(node);
        node->index = singles ++;
        node->is_single = true;
        node->weight = 1;
        node->size = 1;
    }
}

Taxa::Taxa(const Taxa &source) {
    singles = source.singles;
    this->weighting = source.weighting;
    for (auto element : source.label2node) {
        Node *new_node = new Node(element.first);
        label2node[element.first] = new_node;
        new_node->index = element.second->index;
        new_node->is_single = element.second->is_single;
        new_node->weight = element.second->weight;
        new_node->size = element.second->size;
    }
    for (auto element : source.label2node) {
        Node *new_node = label2node[element.first];
        if (element.second->parent == NULL)
            new_node->parent = NULL;
        else 
            new_node->parent = label2node[element.second->parent->label];
    }
    for (Node *root : source.roots) {
        roots.push_back(label2node[root->label]);
    }
    for (Node *leaf : source.leaves) {
        leaves.push_back(label2node[leaf->label]);
    }
}

Taxa::~Taxa() {
    for (auto element : label2node) {
        delete element.second;
    }
}

std::string Taxa::to_string() {
    std::string s = "";
    for (Node *root : roots) 
        s += root->label + " ";
    return s;
}

std::string Taxa::info() {
    return std::to_string(single_size()) + "," + std::to_string(multiple_size());
}

int Taxa::get_weighting() {
    return weighting;
}

void Taxa::update(str<void>::set labels, std::string pseudo) {
    Node *root = new Node(pseudo);
    label2node[pseudo] = root;
    roots.push_back(root);
    root->is_single = false;
    root->weight = 1.0 / labels.size();
    root->size = 0;
    for (std::string label : labels) {
        Node *node = label2node[label];
        node->parent = root;
        node->is_single = false;
        root->size += node->size;
    }
    std::vector<Node *> new_roots = std::vector<Node *>();
    for (Node *root : roots) {
        root->index = -1;
        if (root->parent == NULL && root->is_single) 
            new_roots.push_back(root);
    }
    singles = new_roots.size();
    for (Node *root : roots) {
        root->index = -1;
        if (root->parent == NULL && ! root->is_single) 
            new_roots.push_back(root);
    }
    roots.clear();
    int i = 0;
    for (Node *root : new_roots) {
        roots.push_back(root);
        root->index = i ++;
    }
}

int Taxa::leaf_size() {
    return leaves.size();
}

int Taxa::size() {
    return roots.size();
}

int Taxa::single_size() {
    return singles;
}

int Taxa::multiple_size() {
    return size() - single_size();
}

std::string Taxa::index2leaflabel(int i) {
    return leaves[i]->label;
}

std::string Taxa::index2label(int i) {
    return roots[i]->label;
}

int Taxa::label2index(std::string label) {
    Node *root = get_root(label);
    return root->index;
}

int Taxa::label2key(std::string label) {
    if (is_single(label)) return 0;
    return label2index(label) - singles + 1;
}

std::string Taxa::key2label(int i) {
    return roots[i - 1 + singles]->label;
}

std::string Taxa::get_label(std::string label) {
    Node *root = get_root(label);
    return root->label;
}

double Taxa::get_weight(std::string label) {
    if (weighting == 0) {
        Node *node = label2node[label];
        return get_weight(node);
    }
    else if (weighting == 1) {
        Node *root = get_root(label);
        return 1.0 / root->size;
    }
    else {
        return 1.0;
    }
}

bool Taxa::is_single(std::string label) {
    Node *node = label2node[label];
    return node->is_single;
}

Taxa::Node *Taxa::get_root(std::string label) {
    return get_root(label2node[label]);
}

Taxa::Node *Taxa::get_root(Node *node) {
    if (node->parent == NULL) 
        return node;
    return get_root(node->parent);
}

double Taxa::get_weight(Node *node) {
    if (node->parent == NULL) 
        return node->weight;
    return node->weight * get_weight(node->parent);
}

int Tree::Node::pseudonyms = 0;

Tree::Node::Node(const std::string &name) {
    left = right = parent = NULL;
    size = -1;
    leaves = NULL;
    pairs = NULL;
    s1 = s2 = 0;
    label = name;
}

void Tree::Node::new_states(int size) {
    this->size = size;
    leaves = new double[size + 1];
    for (int i = 0; i < size + 1; i ++) 
        leaves[i] = 0;
    pairs = new double*[size + 1];
    for (int i = 0; i < size + 1; i ++) {
        pairs[i] = new double[size + 1];
        for (int j = 0; j < size + 1; j ++) {
            pairs[i][j] = 0;
        }
    }
}

std::string Tree::Node::pseudonym() {
    return "X" + std::to_string(pseudonyms ++);
}

void Tree::Node::delete_states() {
    if (size >= 0) {
        delete [] leaves;
        for (int i = 0; i < size + 1; i ++) 
            delete [] pairs[i];
        delete [] pairs;
        size = -1;
    }
}

Tree::Node::~Node() {
    delete left;
    delete right;
}

double Tree::Node::get_pairs(double *leaves, double s1, double s2, int x, int y) {
    double t0 = leaves[0], t1 = s1, t2 = s2;
    if (x != 0) {
        t1 -= leaves[x];
        t2 -= leaves[x] * leaves[x];
    }
    if (y != 0) {
        t1 -= leaves[y];
        t2 -= leaves[y] * leaves[y];
    }
    return (t1 * t1 - t2) / 2 + t1 * t0 + t0 * (t0 - 1) / 2;
}

Tree::Tree(std::ifstream &fin, int execution, int weighting, int s0, int s1) {
    std::string newick;
    std::vector<Tree *> input;
    str<void>::set labels;
    srand(s0);
    while (std::getline(fin, newick)) {
        Tree *t = new Tree(newick);
        input.push_back(t);
        t->append_labels(t->root, labels);
    }
    Taxa subset = Taxa(labels, weighting);
    srand(s1);
    if (execution == 0) {
        root = construct_stree(input, subset);
    }
    else if (execution == 1) {
        str<double>::map quartets;
        for (Tree *t : input) t->append_quartets(quartets, subset);
        root = construct_stree_brute(quartets, subset);
    }
    else {
        str<double>::map quartets;
        for (Tree *t : input) t->append_quartets(quartets, subset);
        root = construct_stree_check(quartets, input, subset);
    }
    for (Tree *t : input) delete t;
}

Tree::Tree(const std::string &newick) {
    root = build_tree(newick);
}

std::string Tree::to_string() {
    return display_tree(root);
}

Tree::~Tree() {
    delete root;
}

void Tree::append_labels(Node *root, str<void>::set &labels) {
    if (root->left == NULL) {
        if (labels.find(root->label) == labels.end()) 
            labels.insert(root->label);
    }
    else {
        append_labels(root->left, labels);
        append_labels(root->right, labels);
    }
}

double ***Tree::build_graph(Taxa &subset) {
    int s = subset.single_size(), m = subset.multiple_size();
    build_states(root, subset);
    for (int i = 0; i <= m; i ++) 
        single_pairs(root, 0, i);
    for (int i = 1; i <= m; i ++) 
        for (int j = 0; j <= m; j ++) 
            if (i != j) multiple_pairs(root, i, j);
    double ***mat = Matrix::new_mat<double>(subset.size());
    build_mat(root, subset, mat);
    if (subset.get_weighting() == 2) {
        double *c = new double[m + 1];
        for (int i = 0; i <= m; i ++) 
            c[i] = root->leaves[i];
        c[0] -= 2;
        double sum = Node::get_pairs(c, root->s1, root->s2, 0, 0);
        for (int i = 0; i < s; i ++) {
            for (int j = 0; j < s; j ++) {
                if (i == j) continue;
                mat[i][j][0] = sum - mat[i][j][1];
            }
        }
        for (int i = 0; i <= m; i ++) 
            c[i] = root->leaves[i];
        c[0] -= 1;
        for (int i = 0; i < s; i ++) {
            for (int j = 0; j < m; j ++) {
                double sum = Node::get_pairs(c, root->s1, root->s2, j + 1, 0) * root->leaves[j + 1];
                mat[s + j][i][0] = mat[i][s + j][0] = sum - mat[i][s + j][1];
            }
        }
        for (int i = 0; i <= m; i ++) 
            c[i] = root->leaves[i];
        for (int i = 0; i < m; i ++) {
            for (int j = 0; j < m; j ++) {
                if (i == j) continue;
                double sum = Node::get_pairs(c, root->s1, root->s2, i + 1, j + 1) * root->leaves[i + 1] * root->leaves[j + 1];
                mat[s + j][s + i][0] = mat[s + i][s + j][0] = sum - mat[s + i][s + j][1];
            }
        }
        delete [] c;
    }
    else {
        double *c = new double[m + 1];
        for (int i = 0; i <= m; i ++) 
            c[i] = i == 0 ? root->leaves[i] - 2 : root->leaves[i];
        double sum = Node::get_pairs(c, root->s1, root->s2, 0, 0);
        delete [] c;
        for (int i = 0; i < subset.size(); i ++) {
            for (int j = 0; j < subset.size(); j ++) {
                if (i == j) continue;
                mat[i][j][0] = sum - mat[i][j][1];
            }
        }
    }
    clear_states(root);
    return mat;
}

double Tree::pairs(Node *subtree, int x, int y, bool comple) {
    if (! comple) return Node::get_pairs(subtree->leaves, subtree->s1, subtree->s2, x, y);
    double *c = new double[subtree->size + 1], s1 = 0, s2 = 0;
    for (int i = 0; i <= subtree->size; i ++) {
        c[i] = root->leaves[i] - subtree->leaves[i];
        if (i != 0) {
            s1 += c[i];
            s2 += c[i] * c[i];
        }
    }
    double ret = Node::get_pairs(c, s1, s2, x, y);
    delete [] c;
    return ret;
}

Tree::Node *Tree::build_tree(const std::string &newick) {
    if (newick.length() == 0 || newick.at(0) != '(') {
        std::string delimiter = ":";
        std::string label = newick.substr(0, newick.find(delimiter));
        Node *root = new Node(label);
        label2node[label] = root;
        return root;
    }
    else {
        std::vector<Node *> subtrees;
        int k = 1;
        for (int i = 0, j = 0; i < newick.length(); i ++) {
            if (newick.at(i) == '(') j ++;
            if (newick.at(i) == ')') j --;
            if (newick.at(i) == ',' && j == 1) {
                subtrees.push_back(build_tree(newick.substr(k, i - k)));
                k = i + 1;
            }
        }
        int i = newick.length() - 1;
        while (newick.at(i) != ')') i --;
        subtrees.push_back(build_tree(newick.substr(k, i - k)));
        while (subtrees.size() > 1) {
            int i = rand() % subtrees.size(), j = i; 
            while (j == i) j = rand() % subtrees.size();
            Node *root = new Node(Node::pseudonym());
            root->left = subtrees[i]; root->right = subtrees[j];
            root->left->parent = root->right->parent = root;
            subtrees.erase(subtrees.begin() + i);
            subtrees.erase(subtrees.begin() + (j > i ? j - 1 : j));
            subtrees.push_back(root);
        }
        return subtrees[0];
    }
}

std::string Tree::display_tree(Node *root) {
    if (root->left == NULL) 
        return root->label;
    return "(" + display_tree(root->left) + "," + display_tree(root->right) + ")";
}

void Tree::clear_states(Node *root) {
    root->delete_states();
    if (root->left != NULL) {
        clear_states(root->left);
        clear_states(root->right);
    }
}

void Tree::build_states(Node *root, Taxa &subset) {
    root->new_states(subset.multiple_size());
    if (root->left == NULL) {
        int key = subset.label2key(subset.get_label(root->label));
        root->leaves[key] = subset.get_weight(root->label);
    }
    else {
        build_states(root->left, subset);
        build_states(root->right, subset);
        for (int i = 0; i < root->size + 1; i ++) 
            root->leaves[i] = root->left->leaves[i] + root->right->leaves[i];
    }
    root->s1 = root->s2 = 0;
    for (int i = 1; i < root->size + 1; i ++) {
        root->s1 += root->leaves[i];
        root->s2 += root->leaves[i] * root->leaves[i];
    }
}

void Tree::build_depth(Node *root, int depth) {
    root->new_states(0);
    root->pairs[0][0] = depth;
    if (root->left != NULL) {
        build_depth(root->left, depth + 1);
        build_depth(root->right, depth + 1);
    }
}

void Tree::single_pairs(Node *root, double s, int x) {
    root->pairs[0][x] = s;
    if (root->left != NULL) {
        single_pairs(root->left, s + pairs(root->right, x, 0, false), x);
        single_pairs(root->right, s + pairs(root->left, x, 0, false), x);
    }
}

double Tree::multiple_pairs(Node *root, int x, int y) {
    if (root->left == NULL) return 0;
    double t0 = pairs(root->left, x, y, false), t1 = pairs(root->right, x, y, false);
    double s0 = t0 == 0 ? 0 : multiple_pairs(root->left, x, y);
    double s1 = t1 == 0 ? 0 : multiple_pairs(root->right, x, y);
    double s = s0 + s1 + root->left->leaves[x] * t1 + root->right->leaves[x] * t0;
    root->pairs[x][y] = s;
    return s;
}

str<void>::set Tree::build_mat(Node *root, Taxa &subset, double ***mat) {
    if (root->left == NULL) {
        str<void>::set subtree;
        std::string label = subset.get_label(root->label);
        subtree.insert(label);
        return subtree;
    }
    else {
        str<void>::set left = build_mat(root->left, subset, mat);
        str<void>::set right = build_mat(root->right, subset, mat);
        str<void>::set subtree;
        for (auto i = left.begin(); i != left.end(); i ++) {
            if (subtree.find(*i) == subtree.end()) 
                subtree.insert(*i);
        }
        for (auto j = right.begin(); j != right.end(); j ++) {
            if (subtree.find(*j) == subtree.end()) 
                subtree.insert(*j);
        }
        for (auto i = left.begin(); i != left.end(); i ++) {
            for (auto j = right.begin(); j != right.end(); j ++) {
                if (*i == *j) continue;
                int x = subset.label2key(*i), y = subset.label2key(*j);
                int i_ = subset.label2index(*i), j_ = subset.label2index(*j);
                if (x == 0) {
                    if (y == 0) {
                        double s = 0;
                        s += label2node[*i]->pairs[0][0] + label2node[*j]->pairs[0][0];
                        s -= root->left->pairs[0][0] + root->right->pairs[0][0];
                        s += pairs(root, 0, 0, true);
                        mat[i_][j_][1] = mat[j_][i_][1] = s;
                    }
                    else {
                        double s = root->right->pairs[y][0], t = root->right->leaves[y];
                        s += t * (label2node[*i]->pairs[0][y] - root->left->pairs[0][y]);
                        s += t * pairs(root, y, 0, true);
                        mat[i_][j_][1] += s;
                        mat[j_][i_][1] += s;
                    }
                }
                else {
                    if (y == 0) {
                        double s = root->left->pairs[x][0], t = root->left->leaves[x];
                        s += t * (label2node[*j]->pairs[0][x] - root->right->pairs[0][x]);
                        s += t * pairs(root, x, 0, true);
                        mat[i_][j_][1] += s;
                        mat[j_][i_][1] += s;
                    }
                    else {
                        double s = 0, t0 = root->left->leaves[x], t1 = root->right->leaves[y];
                        s += t1 * root->left->pairs[x][y] + t0 * root->right->pairs[y][x];
                        s += t1 * t0 * pairs(root, x, y, true);
                        mat[i_][j_][1] += s;
                        mat[j_][i_][1] += s;
                    }
                }
            }
        }
        return subtree;
    }
}

Tree::Node *Tree::reroot(Node *root, str<void>::set &visited) {
    std::vector<Node *> child;
    if (root->parent != NULL && visited.find(root->parent->label) == visited.end()) {
        visited.insert(root->parent->label);
        child.push_back(reroot(root->parent, visited));
    }
    if (root->left != NULL && visited.find(root->left->label) == visited.end()) {
        visited.insert(root->left->label);
        child.push_back(reroot(root->left, visited));
    }
    if (root->right != NULL && visited.find(root->right->label) == visited.end()) {
        visited.insert(root->right->label);
        child.push_back(reroot(root->right, visited));
    }
    if (child.size() == 2) {
        Node *new_root = new Node(Node::pseudonym());
        visited.insert(new_root->label);
        new_root->left = child[0];
        new_root->right = child[1];
        new_root->left->parent = new_root->right->parent = new_root;
        return new_root;
    }
    else if (child.size() == 1) {
        return child[0];
    }
    else {
        Node *new_root = new Node(root->label);
        return new_root;
    }
}

Tree::Node *Tree::pseudo2node(Node *root, const std::string &pseudo) {
    if (root->left == NULL) {
        if (root->label == pseudo) 
            return root;
        return NULL;
    }
    else {
        Node *left = pseudo2node(root->left, pseudo);
        if (left != NULL) return left;
        Node *right = pseudo2node(root->right, pseudo);
        if (right != NULL) return right;
        return NULL;
    }
}

Tree::Node *Tree::reroot_stree(Node *root, const std::string &pseudo) {
    Node *new_root = pseudo2node(root, pseudo);
    str<void>::set visited;
    visited.insert(new_root->label);
    Node *new_tree = reroot(new_root, visited);
    delete root;
    return new_tree;
}

Tree::Node *Tree::construct_stree(std::vector<Tree *> &input, Taxa &subset) {
    stats[0] += 1;
    stats[1] += subset.single_size();
    stats[2] += subset.multiple_size();
    int size = subset.size();
    if (verbose >= 1) std::cout << subset.info() << std::endl;
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.index2label(0));
        }
        else if (size == 2) {
            root = new Node(Node::pseudonym());
            root->left = new Node(subset.index2label(0));
            root->right = new Node(subset.index2label(1));
            root->left->parent = root->right->parent = root;
        }
        else {
            root = new Node(Node::pseudonym());
            root->left = new Node(Node::pseudonym());
            root->left->left = new Node(subset.index2label(0));
            root->left->right = new Node(subset.index2label(1));
            root->left->left->parent = root->left->right->parent = root->left;
            root->right = new Node(subset.index2label(2));
            root->left->parent = root->right->parent = root;
        }
    }
    else {
        Graph *g = new Graph(input, subset);
        str<void>::set As, Bs;
        double weight = g->get_cut(& As, & Bs);
        Taxa Am(subset), Bm(subset);
        std::string pseudo = Tree::Node::pseudonym();
        Bm.update(As, pseudo);
        Am.update(Bs, pseudo);
        root = new Node(Node::pseudonym());
        root->left = reroot_stree(construct_stree(input, Am), pseudo);
        root->right = reroot_stree(construct_stree(input, Bm), pseudo);
        root->left->parent = root->right->parent = root;
        delete g;
    }
    if (verbose >= 2) std::cout << display_tree(root) << std::endl;
    return root;
}

void Tree::append_quartet(str<double>::map &quartets, std::string quartet, double weight) {
    if (quartets.find(quartet) == quartets.end()) 
        quartets[quartet] = 0;
    quartets[quartet] += weight;
}

Tree::Node *Tree::construct_stree_brute(str<double>::map &input, Taxa &subset) {
    stats[0] += 1;
    stats[1] += subset.single_size();
    stats[2] += subset.multiple_size();
    int size = subset.size();
    if (verbose >= 1) std::cout << subset.info() << std::endl;
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.index2label(0));
        }
        else if (size == 2) {
            root = new Node(Node::pseudonym());
            root->left = new Node(subset.index2label(0));
            root->right = new Node(subset.index2label(1));
            root->left->parent = root->right->parent = root;
        }
        else {
            root = new Node(Node::pseudonym());
            root->left = new Node(Node::pseudonym());
            root->left->left = new Node(subset.index2label(0));
            root->left->right = new Node(subset.index2label(1));
            root->left->left->parent = root->left->right->parent = root->left;
            root->right = new Node(subset.index2label(2));
            root->left->parent = root->right->parent = root;
        }
    }
    else {
        Graph *g = new Graph(input, subset);
        str<void>::set As, Bs;
        double weight = g->get_cut(& As, & Bs);
        Taxa Am(subset), Bm(subset);
        std::string pseudo = Tree::Node::pseudonym();
        Bm.update(As, pseudo);
        Am.update(Bs, pseudo);
        str<double>::map Ai, Bi;
        for (auto quartet : input) {
            std::string *labels = Tree::split(quartet.first);
            int j = 0;
            for (int i = 0; i < 4; i ++) 
                if (As.find(subset.get_label(labels[i])) != As.end()) 
                    j ++;
            delete [] labels;
            if (j < 2) append_quartet(Bi, quartet.first, quartet.second);
            if (j > 2) append_quartet(Ai, quartet.first, quartet.second);
        }
        root = new Node(Node::pseudonym());
        root->left = reroot_stree(construct_stree_brute(Ai, Am), pseudo);
        root->right = reroot_stree(construct_stree_brute(Bi, Bm), pseudo);
        root->left->parent = root->right->parent = root;
        delete g;
    }
    if (verbose >= 2) std::cout << display_tree(root) << std::endl;
    return root;
}

Tree::Node *Tree::construct_stree_check(str<double>::map &input, std::vector<Tree *> &input_, Taxa &subset) {
    stats[0] += 1;
    stats[1] += subset.single_size();
    stats[2] += subset.multiple_size();
    int size = subset.size();
    if (verbose >= 1) std::cout << subset.info() << std::endl;
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.index2label(0));
        }
        else if (size == 2) {
            root = new Node(Node::pseudonym());
            root->left = new Node(subset.index2label(0));
            root->right = new Node(subset.index2label(1));
            root->left->parent = root->right->parent = root;
        }
        else {
            root = new Node(Node::pseudonym());
            root->left = new Node(Node::pseudonym());
            root->left->left = new Node(subset.index2label(0));
            root->left->right = new Node(subset.index2label(1));
            root->left->left->parent = root->left->right->parent = root->left;
            root->right = new Node(subset.index2label(2));
            root->left->parent = root->right->parent = root;
        }
    }
    else {
        Graph *g = new Graph(input, subset);
        Graph *g_ = new Graph(input_, subset);
        double diff = g->distance(g_, 0) + g->distance(g_, 1);
        if (diff > 1e-6) std::cout << subset.to_string() << ": " << g->distance(g_, 0) << std::endl;
        str<void>::set As, Bs;
        double weight = g->get_cut(& As, & Bs);
        Taxa Am(subset), Bm(subset);
        std::string pseudo = Tree::Node::pseudonym();
        Bm.update(As, pseudo);
        Am.update(Bs, pseudo);
        str<double>::map Ai, Bi;
        for (auto quartet : input) {
            std::string *labels = Tree::split(quartet.first);
            int j = 0;
            for (int i = 0; i < 4; i ++) 
                if (As.find(subset.get_label(labels[i])) != As.end()) 
                    j ++;
            delete [] labels;
            if (j < 2) append_quartet(Bi, quartet.first, quartet.second);
            if (j > 2) append_quartet(Ai, quartet.first, quartet.second);
        }
        root = new Node(Node::pseudonym());
        root->left = reroot_stree(construct_stree_check(Ai, input_, Am), pseudo);
        root->right = reroot_stree(construct_stree_check(Bi, input_, Bm), pseudo);
        root->left->parent = root->right->parent = root;
        delete g;
    }
    if (verbose >= 2) std::cout << display_tree(root) << std::endl;
    return root;
}

std::string Tree::ordered(const std::string &a, const std::string &b) {
    return a < b ? a + " " + b : b + " " + a;
}

std::string Tree::join(std::string *labels) {
    return ordered(ordered(labels[0], labels[1]), ordered(labels[2], labels[3]));
}

std::string *Tree::split(const std::string &qt) {
    std::string *labels = new std::string[4];
    std::string delimiter = " ", s = qt;
    for (int i = 0, j = 0; i < 4; i ++) {
        if (i == 3) {
            labels[i] = s;
        }
        else {
            j = s.find(delimiter);
            labels[i] = s.substr(0, j);
            s.erase(0, j + delimiter.length());
        }
    }
    return labels;
}

void Tree::append_quartets(str<double>::map &quartets, Taxa &subset) {
    build_depth(root, 0);
    int idx[4], leaves = subset.leaf_size();
    for (idx[0] = 0; idx[0] < leaves; idx[0] ++) 
    for (idx[1] = idx[0] + 1; idx[1] < leaves; idx[1] ++) 
    for (idx[2] = idx[1] + 1; idx[2] < leaves; idx[2] ++) 
    for (idx[3] = idx[2] + 1; idx[3] < leaves; idx[3] ++) {
        Node *nodes[4];
        int deepest = -1, id[4];
        std::string labels[4];
        for (int i = 0; i < 4; i ++) {
            nodes[i] = label2node[subset.index2leaflabel(idx[i])];
            id[i] = -1;
        }
        for (int i = 0; i < 4; i ++) {
            for (int j = i + 1; j < 4; j ++) {
                Node *p = nodes[i], *q = nodes[j];
                while (p->pairs[0][0] > q->pairs[0][0]) p = p->parent;
                while (p->pairs[0][0] < q->pairs[0][0]) q = q->parent;
                while (p != q) { p = p->parent; q = q->parent; }
                if (p->pairs[0][0] > deepest) {
                    deepest = p->pairs[0][0];
                    id[0] = i; id[1] = j;
                }
            }
        }
        for (int i = 0; i < 4; i ++) 
            if (i != id[0] && i != id[1]) 
                id[(id[2] == -1 ? 2 : 3)] = i;
        for (int i = 0; i < 4; i ++) 
            labels[i] = nodes[id[i]]->label;
        append_quartet(quartets, join(labels), 1.0);
    }
    clear_states(root);
}

Graph::Graph(std::vector<Tree*> &input, Taxa &subset) {
    size = subset.size();
    for (int i = 0; i < size; i ++) {
        label2index[subset.index2label(i)] = i;
        labels.push_back(subset.index2label(i));
    }
    graph = Matrix::new_mat<double>(size);
    for (int i = 0; i < input.size(); i ++) {
        Tree *t = input[i];
        double ***mat = t->build_graph(subset);
        for (int j = 0; j < size; j ++) {
            for (int k = 0; k < size; k ++) {
                graph[j][k][0] += mat[j][k][0];
                graph[j][k][1] += mat[j][k][1];
            }
        }
        Matrix::delete_mat<double>(mat, size);
    }
    if (verbose >= 3) Matrix::display_mat<double>(graph, size, 0);
}

Graph::Graph(str<double>::map &input, Taxa &subset) {
    size = subset.size();
    for (int i = 0; i < size; i ++) {
        label2index[subset.index2label(i)] = i;
        labels.push_back(subset.index2label(i));
    }
    graph = Matrix::new_mat<double>(size);
    for (auto quartet : input) {
        std::string *labels = Tree::split(quartet.first);
        int a = label2index[subset.get_label(labels[0])],
            b = label2index[subset.get_label(labels[1])],
            c = label2index[subset.get_label(labels[2])],
            d = label2index[subset.get_label(labels[3])];
        double w = quartet.second;
        for (int i = 0; i < 4; i ++) 
            w *= subset.get_weight(labels[i]);
        graph[a][b][1] += w; graph[c][d][1] += w; graph[b][a][1] += w; graph[d][c][1] += w;
        graph[a][c][0] += w; graph[a][d][0] += w; graph[b][c][0] += w; graph[b][d][0] += w;
        graph[c][a][0] += w; graph[d][a][0] += w; graph[c][b][0] += w; graph[d][b][0] += w;
        delete [] labels;
    }
    if (verbose >= 3) Matrix::display_mat<double>(graph, size, 0);
}

double Graph::distance(Graph *g, int k) {
    if (size != g->size) return -1;
    double sum = 0;
    for (int i = 0; i < size; i ++) {
        for (int j = 0; j < size; j ++) {
            sum += abs(graph[i][j][k] - g->graph[i][j][k]);
        }
    }
    return sum / size / size;
}

double Graph::get_cut(str<void>::set *A, str<void>::set *B) {
    double positive_weight = -1.0;
    str<void>::set a, b;
    double lower = 0.0, upper = 6.0;
    while (lower + 0.1 < upper) {
        double alpha = (lower + upper) / 2.0;
        a.clear(); b.clear();
        double weight = sdp_cut(alpha, &a, &b);
        if (weight < 0.001) {
            upper = alpha;
        }
        else {
            lower = alpha;
            positive_weight = weight;
            *A = a;
            *B = b;
        }
    }
    return positive_weight;
}

Graph::~Graph() {
    Matrix::delete_mat<double>(graph, size);
}

double Graph::sdp_cut(double alpha, str<void>::set *A, str<void>::set *B) {
    std::vector<Instance::InstanceTuple> input;
    double sum = 0;
    for (int i = 0; i < size; i ++) {
        for (int j = i + 1; j < size; j ++) {
            double weight = graph[i][j][0] - alpha * graph[i][j][1];
            if (weight < 0) weight = - weight;
            sum += weight;
        }
    }
    double norm = size * (size - 1) / 2 / sum;
    for (int i = 0; i < size; i ++) {
        for (int j = i + 1; j < size; j ++) {
            double weight = (graph[i][j][0] - alpha * graph[i][j][1]) * norm;
            input.push_back(Instance::InstanceTuple(std::make_pair(i + 1, j + 1), weight));
        }
    }
    MaxCutInstance instance(input, size);
    Burer2002 heuristic(instance, -1, false, NULL);
    MaxCutSimpleSolution solution = heuristic.get_best_solution();
    std::vector<int> cut = solution.get_assignments();
    for (int i = 0; i < cut.size(); i ++) {
        if (cut[i] < 0) 
            A->insert(labels[i]);
        else 
            B->insert(labels[i]);
    }
    return solution.get_weight();
}

int main(int argc, char** argv) {
    std::ifstream fin;
    std::ofstream fout;
    std::string item;
    int execution = 0, weighting = 0;
    for (int i = 0; i < argc; i ++) {
        std::string opt(argv[i]);
        if (opt == "-h") item = i == argc - 1 ? "-h" : argv[++ i];
        if (opt == "-i" && i < argc - 1) fin.open(argv[++ i]);
        if (opt == "-o" && i < argc - 1) fout.open(argv[++ i]);
        if (opt == "-v" && i < argc - 1) verbose = std::stoi(argv[++ i]);
        if (opt == "-x" && i < argc - 1) execution = std::stoi(argv[++ i]);
        if (opt == "-w" && i < argc - 1) weighting = std::stoi(argv[++ i]);
    }
    if (item != "") {
        if (item == "-h") {
            std::cout << "use '-h [option]' for more information" << std::endl;
            std::cout << "options: -h -i -o -v -x -w" << std::endl;
            std::cout << "e.g., ./tQMC -h -i" << std::endl;
        }
        if (item == "-i") {
            std::cout << "use '-i [filename]' to assign input file" << std::endl;
            std::cout << "e.g., ./tQMC -i input.tre -o output.tre" << std::endl;
        }
        if (item == "-o") {
            std::cout << "use '-o [filename]' to assign output file" << std::endl;
            std::cout << "e.g., ./tQMC -i input.tre -o output.tre" << std::endl;
        }
        if (item == "-v") {
            std::cout << "use '-v [option]' to show intermediate results" << std::endl;
            std::cout << "options: 0 = none (default), 1 = subproblems, 2 = subproblems and subtrees, 3 = subproblems, subtrees and graphs" << std::endl;
            std::cout << "e.g., ./tQMC -i input.tre -o output.tre -v 1" << std::endl;
        }
        if (item == "-x") {
            std::cout << "use '-x [option]' to set the execution mode" << std::endl;
            std::cout << "options: 0 = fast (default), 1 = naive, 2 = both and compare the results (test)" << std::endl;
            std::cout << "e.g., ./tQMC -i input.tre -o output.tre -x 1" << std::endl;
        }
        if (item == "-w") {
            std::cout << "use '-w [option]' to set the weighting approach" << std::endl;
            std::cout << "options: 0 = normalization and weighting (default), 1 = only normalization, 2 = none" << std::endl;
            std::cout << "e.g., ./tQMC -i input.tre -o output.tre -w 1" << std::endl;
        }
    }
    else {
        if (! fin.is_open()) {
            std::cout << "input file error" << std::endl;
            std::cout << "try ./tQMC -h for help information" << std::endl;
        }
        else if (! fout.is_open()) {
            std::cout << "output file error" << std::endl;
            std::cout << "try ./tQMC -h for help information" << std::endl;
        }
        else {
            stats[0] = stats[1] = stats[2] = 0;
            auto start = std::chrono::high_resolution_clock::now();
            Tree *t = new Tree(fin, execution, weighting, 1, 1);
            fout << t->to_string() << ';' << std::endl;
            delete t;
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "execution time: " << duration.count() << "ms" << std::endl;
            std::cout << "subproblems: " << stats[0] << std::endl;
            std::cout << "input taxa: " << stats[1] << std::endl;
            std::cout << "artificial taxa: " << stats[2] << std::endl;
        }
    }
    return 0;
}
