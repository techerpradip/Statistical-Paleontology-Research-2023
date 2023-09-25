#include <iostream>
#include <unordered_set>
#include <vector>
#include <math.h>
#include <iomanip>
#include <chrono>
#include "pthread.h"
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;


/**
 * @brief A node structure used to represent both internal nodes and species from a given phylogenetic tree.
 * 
 *        id: Unique ID of a given node in the tree.
 *        rate: Incoming numnber of changes to a given node in the tree.
 *        left: Pointer to the left children node of a given node.
 *        right: Pointer to the right children node of a given node.
 *        parent: Pointer to the parent node of a given node.
 * 
 */
struct Node {
    int id;
    int rate;
    string name;
    Node* left = nullptr;
    Node* right = nullptr;
    Node* parent = nullptr;

    // Constructor of the Node structure.
    Node(int id, int rate, string name): id(id), rate(rate), name(name) {}
};


/**
 * @brief A ThreadData structure used to store information needed for a thread to run in parallel while running bruteforce algorithm.
 * 
 *        tree_nodes: Reference to a vector of all valid nodes of the tree.
 *        final_result: Reference to a vector of vectors where a vector at a given index i is the best (i+1) partition of the given tree.
 *        partition_count: The size of partition for specific thread to work.
 * 
 */
struct ThreadData {
    vector<Node*>& tree_nodes;
    vector<vector<Node*>>& final_result;
    int partition_count;

    // Constructor of the ThreadData structure.
    ThreadData(vector<Node*>& tree_nodes, vector<vector<Node*>>& final_result, int partition_count): tree_nodes(tree_nodes), final_result(final_result), partition_count(partition_count) {}
};


/**
 * @brief Global variables used while running the bruteforce algorithm in parallel.
 * 
 *        global_min: Global mininum AIC/BIC for a given run of bruteforce algorithm.
 *        global_max: Vector of Nodes with the minimum AIC/BIC for a given run of bruteforce algorithm.
 *        my_mutex: Parallel threads mutex to synchronize threads while running in parallel.
 * 
 */
double global_min;
vector<Node*> global_data;
pthread_mutex_t my_mutex;


/**
 * @brief A helper functin for parse_newick() function to split part of the string into name and rate.
 * 
 * @param current_string 
 * @param name 
 * @param rate 
 * @param separator 
 */
void split_string(string& current_string, string& name, int& rate, char separator = ':') {
    int index = 0;
    string tree_name = "";
    string tree_rate = "";
    bool separatorFound = false;

    while (current_string[index] != '\0') {
        if (current_string[index] != separator) {
            if (!separatorFound) {
                tree_name += current_string[index]; 
            } else {
                tree_rate += current_string[index];
            } 
        } else {
            separatorFound = true;
        }
        index++;
    }

    if (tree_name == "") {
        name = "Internal Node ";
    } else {
        name = tree_name;
    }

    if (tree_rate == "") {
        rate = 0;
    } else {
        rate = stoi(tree_rate);
    }
}


/**
 * @brief A parsing function that parses a given phylogenetic tree in Newick Format and creates a Binary Tree represeting the phylogenetic tree.
 * 
 * @param newick_string Phylogenetic tree in newick format to be parsed.
 * @param node_array Reference to a vector of Node pointers in the caller function to store pointers to all the nodes of the tree.
 * @return Pointer to the root of the binary tree created.
 */
Node* parse_newick(string& newick_string, vector<Node*>& node_array) {
    int size = newick_string.length();
    string current_string = "";
    Node* new_node = nullptr;
    Node* root = nullptr;
    Node* current_parent = nullptr;
    int current_id = 1;
    int rate;
    string name;

    string newick_formatted = "";
    for (int index = 0; index < size; index++) {
        if (newick_string[index] != ' ') {
            newick_formatted += newick_string[index];
        }
    }
    size = newick_formatted.length();

    for (int index = size - 2; index >= 0; index--) {
        if (newick_formatted[index] == ')') {
            split_string(current_string, name, rate);
            new_node = new Node(current_id, rate, name);
            node_array.push_back(new_node);

            if (current_id == 1) {
                root = new_node;
            } else {
                new_node->parent = current_parent;
                if (current_parent->left == nullptr) {
                    current_parent->left = new_node;
                } else {
                    current_parent->right = new_node;
                }
            }

            current_parent = new_node;
            current_string = "";
            current_id++;
        } else if (newick_formatted[index] == '(') {
            if (current_string != "") {
                split_string(current_string, name, rate, ':');
                new_node = new Node(current_id, rate, name);
                node_array.push_back(new_node);

                new_node->parent = current_parent;
                if (current_parent->left == nullptr) {
                    current_parent->left = new_node;
                } else {
                    current_parent->right = new_node;
                }
            }
            Node* temp = current_parent->parent;
            current_parent = temp;
            current_string = "";
            current_id++;
        } else if (newick_formatted[index] == ',') {
            if (newick_formatted[index+1] != '(') {
                split_string(current_string, name, rate, ':');
                new_node = new Node(current_id, rate, name);
                node_array.push_back(new_node);

                new_node->parent = current_parent;
                if (current_parent->left == nullptr) {
                    current_parent->left = new_node;
                } else {
                    current_parent->right = new_node;
                }

                current_string = "";
                current_id++;
            }
        } else {
            current_string = newick_formatted[index] + current_string;
        }
    }

    return root;
}


/**
 * @brief 
 * 
 * @param node 
 * @param partition_set 
 * @param i 
 * @return int 
 */
int sum_of_rates(Node* node, unordered_set<Node*>& partition_set, int i) {
    if (node == nullptr) {
        return 0;
    }
    if (partition_set.find(node) != partition_set.end() && i == 1) {
        return 0;
    }
    if (node->id == 1) {
        return sum_of_rates(node->left, partition_set, 1) + sum_of_rates(node->right, partition_set, 1);
    }
    return node->rate + sum_of_rates(node->left, partition_set, 1) + sum_of_rates(node->right, partition_set, 1);
}


/**
 * @brief 
 * 
 * @param node 
 * @param partition_set 
 * @param i 
 * @return int 
 */
int edge_count(Node* node, unordered_set<Node*>& partition_set, int i) {
    if (node == nullptr) {
        return 0;
    }
    if (partition_set.find(node) != partition_set.end() && i == 1) {
        return 0;
    }
    if (node->id == 1) {
        return edge_count(node->left, partition_set, 1) + edge_count(node->right, partition_set, 1);
    }
    return 1 + edge_count(node->left, partition_set, 1) + edge_count(node->right, partition_set, 1);
}


/**
 * @brief 
 * 
 * @param x 
 * @param lambdahat 
 * @return double 
 */
double log_dpois(int x, double lambdahat) {
    return  (log(pow(lambdahat, x) * exp(-lambdahat) / (tgamma(x+1))));
}

double calculate_likelihood(Node* node, double average, unordered_set<Node*>& partition_set, int i) {
    if (node == nullptr) {
        return 0;
    }
    if (partition_set.find(node) != partition_set.end() && i == 1) {
        return 0;
    }
    if (node->id == 1) {
        return calculate_likelihood(node->left, average, partition_set, 1) + calculate_likelihood(node->right, average, partition_set, 1);
    }
    return log_dpois(node->rate, average) + calculate_likelihood(node->left, average, partition_set, 1) + calculate_likelihood(node->right, average, partition_set, 1);
}


/**
 * @brief 
 * 
 * @param combination 
 * @return double 
 */
double calculate_AIC(vector<Node*>& combination) {
    Node* node;
    int sum; int edges; double average; double aic_value;
    double likelihood = 0;
    int k = combination.size();
    unordered_set<Node*> partition_set(combination.begin(), combination.end());

    for (int i = 0; i < k; i++) {
        node = combination[i];
        sum = sum_of_rates(node, partition_set, 0);
        edges = edge_count(node, partition_set, 0);
        if (edges == 0) {
            average = 0;
        } else {
            average = double(sum) / double(edges);
        }
        likelihood += calculate_likelihood(node, average, partition_set, 0);
    }
    aic_value = 2 * (k - likelihood);
    return aic_value;
}


/**
 * @brief 
 * 
 * @param combination 
 * @param nbranches 
 * @return double 
 */
double calculate_BIC(vector<Node*>& combination, int nbranches) {
    Node* node;
    int sum; int edges; double average; double bic_value;
    double likelihood = 0;
    int k = combination.size();
    unordered_set<Node*> partition_set(combination.begin(), combination.end());

    for (int i = 0; i < k; i++) {
        node = combination[i];
        sum = sum_of_rates(node, partition_set, 0);
        edges = edge_count(node, partition_set, 0);
        if (edges == 0) {
            average = 0;
        } else {
            average = double(sum) / double(edges);
        }
        likelihood += calculate_likelihood(node, average, partition_set, 0);
    }
    bic_value = k * log(nbranches) - 2 * likelihood;
    return bic_value;
}


/**
 * @brief 
 * 
 * @param nodes 
 * @return vector<string> 
 */
vector<string> node_ids(vector<Node*>& nodes) {
    vector<string> result;
    for (int i = 0; i < nodes.size(); i++) {
        result.push_back(nodes[i]->name);
    }
    return result;
}


/**
 * @brief 
 * 
 * @param values 
 */
void print_vector(vector<string>& values) {
    int i;
    cout << "[";
    for (i = 0; i < values.size() - 1; i++) {
        cout << values[i] << ", ";
    }
    cout << values[i] << "] ";
}


/**
 * @brief 
 * 
 * @param nums 
 * @param final_result 
 * @param k 
 * @param start 
 * @param end 
 * @param index 
 * @param current 
 * @param min_aic 
 */
void generate_combinations_aic(vector<Node*>& nums, vector<vector<Node*>>& final_result, int k, int start, int end, int index, vector<Node*>& current, double& min_aic) {
    if (index == k) {
        double aic = calculate_AIC(current);
        if (aic < min_aic) {
            min_aic = aic;
            final_result[k - 1] = current;
        }
        return;
    }
    for (int i = start; i <= end && end - i + 1 >= k - index; i++) {
        current[index] = nums[i];
        generate_combinations_aic(nums, final_result, k, i+1, end, index+1, current, min_aic);
    }
}


/**
 * @brief 
 * 
 * @param arg 
 * @return void* 
 */
void* process_tree_parallel_aic(void* arg) {
    ThreadData current_thread = *(ThreadData*) arg;
    int partition_count = current_thread.partition_count;
    vector<Node*>& nodes_array = current_thread.tree_nodes;
    vector<vector<Node*>>& final_result = current_thread.final_result;

    vector<Node*> current(partition_count, nullptr);
    double min_aic = INFINITY;
    current[0] = nodes_array[0];
    generate_combinations_aic(nodes_array, final_result, partition_count, 1, nodes_array.size() - 1, 1, current, min_aic);

    pthread_mutex_lock(&my_mutex);
    if (min_aic < global_min) {
        global_min = min_aic;
        global_data = final_result[partition_count - 1];
    }
    pthread_mutex_unlock(&my_mutex);

    vector<string> min_aic_ids = node_ids(final_result[partition_count - 1]);
    print_vector(min_aic_ids);
    cout << setprecision(20) << min_aic << endl;

    pthread_exit(nullptr);
}


/**
 * @brief 
 * 
 * @param nums 
 * @param final_result 
 * @param k 
 * @param start 
 * @param end 
 * @param index 
 * @param current 
 * @param min_bic 
 */
void generate_combinations_bic(vector<Node*>& nums, vector<vector<Node*>>& final_result, int k, int start, int end, int index, vector<Node*>& current, double& min_bic) {
    if (index == k) {
        double bic = calculate_BIC(current, nums.size() - 1);
        if (bic < min_bic) {
            min_bic = bic;
            final_result[k - 1] = current;
        }
        return;
    }
    for (int i = start; i <= end && end - i + 1 >= k - index; i++) {
        current[index] = nums[i];
        generate_combinations_bic(nums, final_result, k, i+1, end, index+1, current, min_bic);
    }
}


/**
 * @brief 
 * 
 * @param arg 
 * @return void* 
 */
void* process_tree_parallel_bic(void* arg) {
    ThreadData current_thread = *(ThreadData*) arg;
    int partition_count = current_thread.partition_count;
    vector<Node*>& nodes_array = current_thread.tree_nodes;
    vector<vector<Node*>>& final_result = current_thread.final_result;

    vector<Node*> current(partition_count, nullptr);
    double min_bic = INFINITY;
    current[0] = nodes_array[0];
    generate_combinations_bic(nodes_array, final_result, partition_count, 1, nodes_array.size() - 1, 1, current, min_bic);

    pthread_mutex_lock(&my_mutex);
    if (min_bic < global_min) {
        global_min = min_bic;
        global_data = final_result[partition_count - 1];
    }
    pthread_mutex_unlock(&my_mutex);

    vector<string> min_bic_ids = node_ids(final_result[partition_count - 1]);
    print_vector(min_bic_ids);
    cout << setprecision(20) << min_bic << endl;

    pthread_exit(nullptr);
}


/**
 * @brief 
 * 
 * @param node 
 */
void delete_tree_nodes(Node* node) {
    if (node != nullptr) {
        delete_tree_nodes(node->left);
        delete_tree_nodes(node->right);
        delete node;
    }
}


/**
 * @brief 
 * 
 * @param currNode 
 * @param cutNodes 
 * @param bestPair 
 * @return pair<double, Node*> 
 */
pair<double, Node*> forwardPartition_aic(Node* currNode, vector<Node*> cutNodes, pair<double, Node*> bestPair) {
    if (currNode->left != nullptr) {
        bestPair = forwardPartition_aic(currNode->left, cutNodes, bestPair);
    }
    if (currNode->right != nullptr) {
        bestPair = forwardPartition_aic(currNode->right, cutNodes, bestPair);
    }
    bool cutNodeFound = 0;
    double tempAIC;
    for (Node* checkNode : cutNodes) {
        if (currNode->name == checkNode->name) {
            cutNodeFound = 1;
            break;
        }
    }
    if (!cutNodeFound) {
        cutNodes.push_back(currNode);
        tempAIC = calculate_AIC(cutNodes);
        if (tempAIC < bestPair.first) {
            bestPair.first = tempAIC;
            bestPair.second = currNode;
        }
        cutNodes.pop_back();
    }
    return bestPair;
}


/**
 * @brief 
 * 
 * @param currNode 
 * @param cutNodes 
 * @param bestPair 
 * @param nbranches 
 * @return pair<double, Node*> 
 */
pair<double, Node*> forwardPartition_bic(Node* currNode, vector<Node*> cutNodes, pair<double, Node*> bestPair, int nbranches) {
    if (currNode->left != nullptr) {
        bestPair = forwardPartition_bic(currNode->left, cutNodes, bestPair, nbranches);
    }
    if (currNode->right != nullptr) {
        bestPair = forwardPartition_bic(currNode->right, cutNodes, bestPair, nbranches);
    }
    bool cutNodeFound = 0;
    double tempBIC;
    for (Node* checkNode : cutNodes) {
        if (currNode->name == checkNode->name) {
            cutNodeFound = 1;
            break;
        }
    }
    if (!cutNodeFound) {
        cutNodes.push_back(currNode);
        tempBIC = calculate_BIC(cutNodes, nbranches);
        if (tempBIC < bestPair.first) {
            bestPair.first = tempBIC;
            bestPair.second = currNode;
        }
        cutNodes.pop_back();
    }
    return bestPair;
}


/**
 * @brief 
 * 
 * @param node_array 
 * @param initial_cuts 
 * @return pair<double, vector<string>> 
 */
pair<double, vector<string>> forwardRecursive_aic(vector<Node*>& node_array, vector<Node*> initial_cuts = {}) {
    const pair<double, Node*> defaultPair = {INFINITY, new Node(999, 999, "Default")};

    vector<Node*> cuts;
    vector<double> AICs;
    pair<double, Node*> tempPair;
    int numNodes;
    int initial_cuts_size;
    double tempAIC;
    pair<double, Node*> bestCutoff;

    numNodes = node_array.size();
    initial_cuts_size = initial_cuts.size();

    if (initial_cuts_size == 0) {
        cuts = {node_array[0]};
        AICs.push_back(calculate_AIC(cuts));
        bestCutoff = {AICs.front(), cuts.front()};
    } else {
        cuts = initial_cuts;
        bestCutoff = {calculate_AIC(cuts), cuts[initial_cuts_size - 1]};
    }

    for (int i = 1; i < numNodes; i++) {
        tempPair = defaultPair; //becuase nothing has been done yet for this iteration
        tempPair = forwardPartition_aic(node_array[0], cuts, tempPair);
        tempAIC = tempPair.first;
        AICs.push_back(tempAIC);
        cuts.push_back(tempPair.second);

        if (tempAIC < bestCutoff.first) {
            bestCutoff.first = tempAIC;
            bestCutoff.second = tempPair.second;
        }
    }

    for (int j = initial_cuts_size; j < cuts.size(); j++) {
        cout<<j<<" cuts ("<<j+1<<" rates): " << cuts.at(j)->name << ": " << setprecision(20) << AICs.at(j - initial_cuts_size) << endl;
    }

   vector<string> retVec;
    for (Node* cutNode : cuts) {
        retVec.push_back(cutNode->name);
        if (cutNode->name == bestCutoff.second->name) {
            break;
        }
    }

    delete defaultPair.second;

    pair<double, vector<string>> retPair = {bestCutoff.first, retVec};
    return retPair;
}


/**
 * @brief 
 * 
 * @param node_array 
 * @param nbranches 
 * @param initial_cuts 
 * @return pair<double, vector<string>> 
 */
pair<double, vector<string>> forwardRecursive_bic(vector<Node*>& node_array, int nbranches, vector<Node*> initial_cuts = {}) {
    const pair<double, Node*> defaultPair = {INFINITY, new Node(999, 999, "Default")};

    vector<Node*> cuts;
    vector<double> BICs;
    pair<double, Node*> tempPair;
    int numNodes;
    int initial_cuts_size;
    double tempBIC;
    pair<double, Node*> bestCutoff;

    numNodes = node_array.size();
    initial_cuts_size = initial_cuts.size();

    if (initial_cuts_size == 0) {
        cuts = {node_array[0]};
        BICs.push_back(calculate_BIC(cuts, nbranches));
        bestCutoff = {BICs.front(), cuts.front()};
    } else {
        cuts = initial_cuts;
        bestCutoff = {calculate_BIC(cuts, nbranches), cuts[initial_cuts_size - 1]};
    }

    for (int i = 1; i < numNodes; i++) {
        tempPair = defaultPair; //becuase nothing has been done yet for this iteration
        tempPair = forwardPartition_bic(node_array[0], cuts, tempPair, nbranches);
        tempBIC = tempPair.first;
        BICs.push_back(tempBIC);
        cuts.push_back(tempPair.second);

        if (tempBIC < bestCutoff.first) {
            bestCutoff.first = tempBIC;
            bestCutoff.second = tempPair.second;
        }
    }

    for (int j = initial_cuts_size; j < cuts.size(); j++) {
        cout<<j<<" cuts ("<<j+1<<" rates): " << cuts.at(j)->name << ": " << setprecision(20) << BICs.at(j - initial_cuts_size) << endl;
    }

   vector<string> retVec;
    for (Node* cutNode : cuts) {
        retVec.push_back(cutNode->name);
        if (cutNode->name == bestCutoff.second->name) {
            break;
        }
    }
    
    delete defaultPair.second;

    pair<double, vector<string>> retPair = {bestCutoff.first, retVec};
    return retPair;
}


/**
 * @brief 
 * 
 * @param currNode 
 * @param root 
 * @param cutNodes 
 * @param bestPair 
 * @return pair<double, Node*> 
 */
pair<double, Node*> backwardMerge_aic(Node* currNode, Node* root, vector<Node*> cutNodes, pair<double, Node*> bestPair) {
     if (currNode->left != nullptr) {
        bestPair = backwardMerge_aic(currNode->left, root, cutNodes, bestPair);
    }
    if (currNode->right != nullptr) {
        bestPair = backwardMerge_aic(currNode->right, root, cutNodes, bestPair);
    }
    bool cutNodeFound = 0;
    double tempAIC;
    int cutSize = cutNodes.size();
    int index;
    if (currNode->name == root->name) {
        return bestPair;
    }
    for (int i = 0; i < cutSize; i++) {
        if (currNode->name == cutNodes[i]->name) {
            cutNodeFound = 1;
            index = i;
            break;
        }
    }
    if (cutNodeFound) {
        Node* tempNode = cutNodes[index];
        cutNodes.erase(cutNodes.begin() + index);
        tempAIC = calculate_AIC(cutNodes);
        if (tempAIC < bestPair.first) {
            bestPair.first = tempAIC;
            bestPair.second = currNode;
        }
        cutNodes.push_back(tempNode);
    }
    return bestPair;
}


/**
 * @brief 
 * 
 * @param currNode 
 * @param root 
 * @param cutNodes 
 * @param bestPair 
 * @param nbranches 
 * @return pair<double, Node*> 
 */
pair<double, Node*> backwardMerge_bic(Node* currNode, Node* root, vector<Node*> cutNodes, pair<double, Node*> bestPair, int nbranches) {
     if (currNode->left != nullptr) {
        bestPair = backwardMerge_bic(currNode->left, root, cutNodes, bestPair, nbranches);
    }
    if (currNode->right != nullptr) {
        bestPair = backwardMerge_bic(currNode->right, root, cutNodes, bestPair, nbranches);
    }
    bool cutNodeFound = 0;
    double tempBIC;
    int cutSize = cutNodes.size();
    int index;
    if (currNode->name == root->name) {
        return bestPair;
    }
    for (int i = 0; i < cutSize; i++) {
        if (currNode->name == cutNodes[i]->name) {
            cutNodeFound = 1;
            index = i;
            break;
        }
    }
    if (cutNodeFound) {
        Node* tempNode = cutNodes[index];
        cutNodes.erase(cutNodes.begin() + index);
        tempBIC = calculate_BIC(cutNodes, nbranches);
        if (tempBIC < bestPair.first) {
            bestPair.first = tempBIC;
            bestPair.second = currNode;
        }
        cutNodes.push_back(tempNode);
    }
    return bestPair;
}


/**
 * @brief 
 * 
 * @param node_array 
 * @return pair<double, vector<string>> 
 */
pair<double, vector<string>> backwardRecursive_aic(vector<Node*>& node_array) {
    const pair<double, Node*> defaultPair = {INFINITY, new Node(999, 999, "Default")};
    
    vector<Node*> cuts(node_array.begin(), node_array.end());
    int numCuts = cuts.size();
    unordered_set<Node*> retSet(cuts.begin(),cuts.end());
    vector<double> AICs;
    
    AICs.push_back(calculate_AIC(cuts));

    pair<double, Node*> tempPair;
    double tempAIC;
    vector<Node*> uncuts = {new Node(999, 999, "None")};
    pair<double, Node*> bestMerge = {AICs.front(), uncuts.front()};
    for (int i = numCuts - 1; i > 0; i--) {
        tempPair = defaultPair; //becuase nothing has been done yet for this iteration
        tempPair = backwardMerge_aic(node_array[0], node_array[0], cuts, tempPair);
        tempAIC = tempPair.first;
        AICs.push_back(tempAIC);
        uncuts.push_back(tempPair.second);
        int index = 0;
        for (int j = 0; j < cuts.size(); j++) {
            if (tempPair.second->name == cuts[j]->name) {
                index = j;
                break;
            }
        }
        cuts.erase(cuts.begin() + index);

        if (tempAIC < bestMerge.first) {
            bestMerge.first = tempAIC;
            bestMerge.second = tempPair.second;
        }
    }

    int m = 0;
    for (Node* uncutNode: uncuts) {
        cout<<numCuts-m-1<<" cuts ("<<numCuts-m<<" rates): " << uncutNode->name << ": " << setprecision(20) << AICs.at(m) << endl;
        m++;
    }

    //Takes all nodes that should be uncut out of the cut set that will be returned
    for (Node* uncutNode : uncuts) {
        retSet.erase(uncutNode);
        if (uncutNode->name == bestMerge.second->name) {
            break;
        }
    }

    vector<string> retVec;
    for (Node* uncutNode : retSet) {
      retVec.push_back(uncutNode->name);
    }

    delete defaultPair.second;
    delete uncuts[0];

    pair<double, vector<string>> retPair = {bestMerge.first, retVec};
    return retPair;
}


/**
 * @brief 
 * 
 * @param node_array 
 * @param nbranches 
 * @return pair<double, vector<string>> 
 */
pair<double, vector<string>> backwardRecursive_bic(vector<Node*>& node_array, int nbranches) {
    const pair<double, Node*> defaultPair = {INFINITY, new Node(999, 999, "Default")};
    
    vector<Node*> cuts(node_array.begin(), node_array.end());
    int numCuts = cuts.size();
    unordered_set<Node*> retSet(cuts.begin(),cuts.end());
    vector<double> BICs;
    
    BICs.push_back(calculate_BIC(cuts, nbranches));

    pair<double, Node*> tempPair;
    double tempBIC;
    vector<Node*> uncuts = {new Node(999, 999, "None")};
    pair<double, Node*> bestMerge = {BICs.front(), uncuts.front()};
    for (int i = numCuts - 1; i > 0; i--) {
        tempPair = defaultPair; //becuase nothing has been done yet for this iteration
        tempPair = backwardMerge_bic(node_array[0], node_array[0], cuts, tempPair, nbranches);
        tempBIC = tempPair.first;
        BICs.push_back(tempBIC);
        uncuts.push_back(tempPair.second);
        int index = 0;
        for (int j = 0; j < cuts.size(); j++) {
            if (tempPair.second->name == cuts[j]->name) {
                index = j;
                break;
            }
        }
        cuts.erase(cuts.begin() + index);

        if (tempBIC < bestMerge.first) {
            bestMerge.first = tempBIC;
            bestMerge.second = tempPair.second;
        }
    }

    int m = 0;
    for (Node* uncutNode: uncuts) {
        cout<<numCuts-m-1<<" cuts ("<<numCuts-m<<" rates): " << uncutNode->name << ": " << setprecision(20) << BICs.at(m) << endl;
        m++;
    }

    //Takes all nodes that should be uncut out of the cut set that will be returned
    for (Node* uncutNode : uncuts) {
        retSet.erase(uncutNode);
        if (uncutNode->name == bestMerge.second->name) {
            break;
        }
    }

    vector<string> retVec;
    for (Node* uncutNode : retSet) {
      retVec.push_back(uncutNode->name);
    }

    delete defaultPair.second;
    delete uncuts[0];
    
    pair<double, vector<string>> retPair = {bestMerge.first, retVec};
    return retPair;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> forwardAlgorithm_aic(string newick_string) {

    vector<Node*> node_array;
    double time_elapsed;
    
    clock_t start = clock();
    
    Node* root = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);

    pair<double, vector<string>> forwardResult;
    cout << "Using Forward Algorithm for AIC:" << endl;
    cout << "" << endl;
    forwardResult = forwardRecursive_aic(node_array);
    print_vector(forwardResult.second);
    cout << "is the " << forwardResult.second.size() << " partition array with the minimum AIC value: " << setprecision(20) << forwardResult.first << endl;
    cout << "" << endl;

    clock_t end = clock();
    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Forward Algorithm for AIC is: " << time_elapsed << endl;
    cout << "" << endl;

    delete_tree_nodes(root);

    return forwardResult.second;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> forwardAlgorithm_bic(string newick_string) {

    vector<Node*> node_array;
    double time_elapsed;
    int nbranches;
    
    clock_t start = clock();
    
    Node* root = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);
    nbranches = node_array.size() - 1;

    pair<double, vector<string>> forwardResult;
    cout << "Using Forward Algorithm for BIC:" << endl;
    cout << "" << endl;
    forwardResult = forwardRecursive_bic(node_array, nbranches);
    print_vector(forwardResult.second);
    cout << "is the " << forwardResult.second.size() << " partition array with the minimum BIC value: " << setprecision(20) << forwardResult.first << endl;
    cout << "" << endl;

    clock_t end = clock();
    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Forward Algorithm for BIC is: " << time_elapsed << endl;
    cout << "" << endl;

    delete_tree_nodes(root);

    return forwardResult.second;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> backwardAlgorithm_aic(string newick_string) {

    vector<Node*> node_array;
    double time_elapsed;
    
    clock_t start = clock();
    
    Node* root = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);

    pair<double, vector<string>> backwardResult;
    cout << "Using Backward Algorithm for AIC:" << endl;
    cout << "" << endl;
    backwardResult = backwardRecursive_aic(node_array);
    print_vector(backwardResult.second);
    cout << "is the " << backwardResult.second.size() << " partition array with the minimum AIC value: " << setprecision(20) << backwardResult.first << endl;
    cout << "" << endl;

    clock_t end = clock();
    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Backward Algorithm for AIC is: " << time_elapsed << endl;
    cout << "" << endl;

    delete_tree_nodes(root);

    return backwardResult.second;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> backwardAlgorithm_bic(string newick_string) {

    vector<Node*> node_array;
    double time_elapsed;
    int nbranches;
    
    clock_t start = clock();
    
    Node* root = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);
    nbranches = node_array.size() - 1;

    pair<double, vector<string>> backwardResult;
    cout << "Using Backward Algorithm for BIC:" << endl;
    cout << "" << endl;
    backwardResult = backwardRecursive_bic(node_array, nbranches);
    print_vector(backwardResult.second);
    cout << "is the " << backwardResult.second.size() << " partition array with the minimum BIC value: " << setprecision(20) << backwardResult.first << endl;
    cout << "" << endl;

    clock_t end = clock();
    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Backward Algorithm for BIC is: " << time_elapsed << endl;
    cout << "" << endl;

    delete_tree_nodes(root);

    return backwardResult.second;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> bruteForceAlgorithm_aic(string newick_string) {

    vector<Node*> node_array;
    double time_elapsed;
    int size;
    pthread_t* threads;

    global_min = INFINITY;
    my_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    clock_t start = clock();
    
    Node* root  = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);
    size = node_array.size();

    cout << "Using Bruteforce Algorithm for AIC:" << endl;
    cout << "" << endl;

    threads = (pthread_t*) malloc(sizeof(pthread_t) * size);
    
    vector<vector<Node*>> final_result(size, {node_array[0]});

    ThreadData current_thread(node_array, final_result, 0);
    vector<ThreadData> threads_information(size, current_thread);


    for (int i = 1; i <= size; i++) {
        threads_information[i-1].partition_count = i;
        int create_result = pthread_create(&threads[i-1], nullptr, process_tree_parallel_aic, &threads_information[i-1]);
        if (create_result != 0) {
            cerr << "Error creating thread!" << endl;
            exit(-1);
        }
    }

    for (int i = 0; i < size; i++) {
        int join_result = pthread_join(threads[i], nullptr);
        if (join_result != 0) {
            cerr << "Error joining thread: " << threads[i] << endl;
            exit(-1);
        }
    }

    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    vector<string> global_min_ids = node_ids(global_data);
    print_vector(global_min_ids);
    cout << "is the " << global_data.size() << " partition array with the minimum AIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Bruteforce Algorithm for AIC is: " << time_elapsed << endl;
    cout << "" << endl;

    free(threads);
    threads = nullptr;

    delete_tree_nodes(root);

    pthread_mutex_destroy(&my_mutex);
    global_data.clear();

    return global_min_ids;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> bruteForceAlgorithm_bic(string newick_string) {

    vector<Node*> node_array;
    double time_elapsed;
    int size;
    pthread_t* threads;

    global_min = INFINITY;
    my_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    clock_t start = clock();
    
    Node* root  = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);
    size = node_array.size();

    cout << "Using Bruteforce Algorithm for BIC:" << endl;
    cout << "" << endl;

    threads = (pthread_t*) malloc(sizeof(pthread_t) * size);
    
    vector<vector<Node*>> final_result(size, {node_array[0]});

    ThreadData current_thread(node_array, final_result, 0);
    vector<ThreadData> threads_information(size, current_thread);


    for (int i = 1; i <= size; i++) {
        threads_information[i-1].partition_count = i;
        int create_result = pthread_create(&threads[i-1], nullptr, process_tree_parallel_bic, &threads_information[i-1]);
        if (create_result != 0) {
            cerr << "Error creating thread!" << endl;
            exit(-1);
        }
    }

    for (int i = 0; i < size; i++) {
        int join_result = pthread_join(threads[i], nullptr);
        if (join_result != 0) {
            cerr << "Error joining thread: " << threads[i] << endl;
            exit(-1);
        }
    }

    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    vector<string> global_min_ids = node_ids(global_data);
    print_vector(global_min_ids);
    cout << "is the " << global_data.size() << " partition array with the minimum BIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Bruteforce Algorithm for BIC is: " << time_elapsed << endl;
    cout << "" << endl;

    free(threads);
    threads = nullptr;

    delete_tree_nodes(root);

    pthread_mutex_destroy(&my_mutex);
    global_data.clear();

    return global_min_ids;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @param brute_count 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> forwardParameterAlgorithm_aic(string newick_string, int brute_count) {

    vector<Node*> node_array;
    double time_elapsed;
    pthread_t* threads;

    global_min = INFINITY;
    my_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    clock_t start = clock();
    
    Node* root  = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);

    cout << "Using Forward with Parameter Algorithm for AIC:" << endl;
    cout << "" << endl;

    threads = (pthread_t*) malloc(sizeof(pthread_t) * brute_count);
    
    vector<vector<Node*>> final_result(brute_count, {node_array[0]});

    ThreadData current_thread(node_array, final_result, 0);
    vector<ThreadData> threads_information(brute_count, current_thread);


    for (int i = 1; i <= brute_count; i++) {
        threads_information[i-1].partition_count = i;
        int create_result = pthread_create(&threads[i-1], nullptr, process_tree_parallel_aic, &threads_information[i-1]);
        if (create_result != 0) {
            cerr << "Error creating thread!" << endl;
            exit(-1);
        }
    }

    for (int i = 0; i < brute_count; i++) {
        int join_result = pthread_join(threads[i], nullptr);
        if (join_result != 0) {
            cerr << "Error joining thread: " << threads[i] << endl;
            exit(-1);
        }
    }

    pair<double, vector<string>> forwardRes = forwardRecursive_aic(node_array, global_data);
    vector<string> global_min_ids;
    if (global_min < forwardRes.first) {
        global_min_ids = node_ids(global_data);
    } else {
        global_min_ids = forwardRes.second;
        global_min = forwardRes.first;
    }

    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    print_vector(global_min_ids);
    cout << "is the " << global_min_ids.size() << " partition array with the minimum AIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Forward with Parameter Algorithm for AIC is: " << time_elapsed << endl;
    cout << "" << endl;

    free(threads);
    threads = nullptr;

    delete_tree_nodes(root);

    pthread_mutex_destroy(&my_mutex);
    global_data.clear();

    return global_min_ids;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @param brute_count 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> forwardParameterAlgorithm_bic(string newick_string, int brute_count) {
    vector<Node*> node_array;
    double time_elapsed;
    int size;
    int nbranches;
    pthread_t* threads;

    global_min = INFINITY;
    my_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    clock_t start = clock();
    
    Node* root  = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);
    size = node_array.size();
    nbranches = size - 1;

    cout << "Using Forward with Parameter Algorithm for BIC:" << endl;
    cout << "" << endl;

    threads = (pthread_t*) malloc(sizeof(pthread_t) * brute_count);
    
    vector<vector<Node*>> final_result(brute_count, {node_array[0]});

    ThreadData current_thread(node_array, final_result, 0);
    vector<ThreadData> threads_information(brute_count, current_thread);


    for (int i = 1; i <= brute_count; i++) {
        threads_information[i-1].partition_count = i;
        int create_result = pthread_create(&threads[i-1], nullptr, process_tree_parallel_bic, &threads_information[i-1]);
        if (create_result != 0) {
            cerr << "Error creating thread!" << endl;
            exit(-1);
        }
    }

    for (int i = 0; i < brute_count; i++) {
        int join_result = pthread_join(threads[i], nullptr);
        if (join_result != 0) {
            cerr << "Error joining thread: " << threads[i] << endl;
            exit(-1);
        }
    }

    pair<double, vector<string>> forwardRes = forwardRecursive_bic(node_array, nbranches, global_data);
    vector<string> global_min_ids;
    if (global_min < forwardRes.first) {
        global_min_ids = node_ids(global_data);
    } else {
        global_min_ids = forwardRes.second;
        global_min = forwardRes.first;
    }

    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    print_vector(global_min_ids);
    cout << "is the " << global_min_ids.size() << " partition array with the minimum BIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Forward with Parameter Algorithm for BIC is: " << time_elapsed << endl;
    cout << "" << endl;

    free(threads);
    threads = nullptr;

    delete_tree_nodes(root);

    pthread_mutex_destroy(&my_mutex);
    global_data.clear();

    return global_min_ids;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @param brute_count 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> backwardParameterAlgorithm_aic(string newick_string, int brute_count) {

    vector<Node*> node_array;
    double time_elapsed;
    int size;
    pthread_t* threads;

    global_min = INFINITY;
    my_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    clock_t start = clock();
    
    Node* root  = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);
    size = node_array.size();

    cout << "Using Backward with Parameter Algorithm for AIC:" << endl;
    cout << "" << endl;

    threads = (pthread_t*) malloc(sizeof(pthread_t) * size);
    
    vector<vector<Node*>> final_result(size, {node_array[0]});

    ThreadData current_thread(node_array, final_result, 0);
    vector<ThreadData> threads_information(size, current_thread);


    for (int i = size; i > size - brute_count; i--) {
        threads_information[i-1].partition_count = i;
        int create_result = pthread_create(&threads[i-1], nullptr, process_tree_parallel_aic, &threads_information[i-1]);
        if (create_result != 0) {
            cerr << "Error creating thread!" << endl;
            exit(-1);
        }
    }

    for (int i = size; i > size - brute_count; i--) {
        int join_result = pthread_join(threads[i - 1], nullptr);
        if (join_result != 0) {
            cerr << "Error joining thread: " << threads[i - 1] << endl;
            exit(-1);
        }
    }

    pair<double, vector<string>> backRes = backwardRecursive_aic(global_data);
    vector<string> global_min_ids;
    if (global_min < backRes.first) {
        global_min_ids = node_ids(global_data);
    } else {
        global_min_ids = backRes.second;
        global_min = backRes.first;
    }

    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    print_vector(global_min_ids);
    cout << "is the " << global_min_ids.size() << " partition array with the minimum AIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Backward with Parameter Algorithm for AIC is: " << time_elapsed << endl;
    cout << "" << endl;

    free(threads);
    threads = nullptr;

    delete_tree_nodes(root);

    pthread_mutex_destroy(&my_mutex);
    global_data.clear();

    return global_min_ids;
}


/**
 * @brief 
 * 
 * @param newick_string 
 * @param brute_count 
 * @return vector<string> 
 */
// [[Rcpp::export]]
vector<string> backwardParameterAlgorithm_bic(string newick_string, int brute_count) {
    vector<Node*> node_array;
    double time_elapsed;
    int size;
    int nbranches;
    pthread_t* threads;

    global_min = INFINITY;
    my_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    clock_t start = clock();
    
    Node* root  = parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);
    size = node_array.size();
    nbranches = size - 1;

    cout << "Using Backward with Parameter Algorithm for BIC:" << endl;
    cout << "" << endl;

    threads = (pthread_t*) malloc(sizeof(pthread_t) * size);
    
    vector<vector<Node*>> final_result(size, {node_array[0]});

    ThreadData current_thread(node_array, final_result, 0);
    vector<ThreadData> threads_information(size, current_thread);


    for (int i = size; i > size - brute_count; i--) {
        threads_information[i-1].partition_count = i;
        int create_result = pthread_create(&threads[i-1], nullptr, process_tree_parallel_bic, &threads_information[i-1]);
        if (create_result != 0) {
            cerr << "Error creating thread!" << endl;
            exit(-1);
        }
    }

    for (int i = size; i > size - brute_count; i--) {
        int join_result = pthread_join(threads[i - 1], nullptr);
        if (join_result != 0) {
            cerr << "Error joining thread: " << threads[i - 1] << endl;
            exit(-1);
        }
    }

    pair<double, vector<string>> backRes = backwardRecursive_bic(global_data, nbranches);
    vector<string> global_min_ids;
    if (global_min < backRes.first) {
        global_min_ids = node_ids(global_data);
    } else {
        global_min_ids = backRes.second;
        global_min = backRes.first;
    }

    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    print_vector(global_min_ids);
    cout << "is the " << global_min_ids.size() << " partition array with the minimum BIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running Backward with Parameter Algorithm for BIC is: " << time_elapsed << endl;
    cout << "" << endl;

    free(threads);
    threads = nullptr;

    delete_tree_nodes(root);

    pthread_mutex_destroy(&my_mutex);
    global_data.clear();

    return global_min_ids;
}
