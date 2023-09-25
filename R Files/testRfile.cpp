#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <iomanip>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

struct Node {
    int id;
    int rate;
    string name;
    Node* left = nullptr;
    Node* right = nullptr;
    Node* parent = nullptr;
    Node(int id, int rate, string name): id(id), rate(rate), name(name) {}
};


void split_string(string& current_string, string& name, int& rate, char separator) {
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


void parse_newick(string& newick_string, vector<Node*>& node_array) {
    int size = newick_string.length();
    string current_string = "";
    Node* new_node = nullptr;
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
            split_string(current_string, name, rate, ':');
            new_node = new Node(current_id, rate, name);
            node_array.push_back(new_node);

            if (current_id != 1) {
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
}

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

void generate_combinations(vector<Node*>& nums, int k, int start, vector<Node*>& current, vector<vector<Node*>>& result) {
    if (k == 0) {
        result.push_back(current);
        return;
    }

    for (int i = start; i < nums.size(); i++) {
        current.push_back(nums[i]);
        generate_combinations(nums, k - 1, i + 1, current, result);
        current.pop_back();
    }
}

vector<vector<Node*>> combinations(vector<Node*>& nums, int k) {
    vector<vector<Node*>> result;
    vector<Node*> current;
    current.push_back(nums[0]);
    generate_combinations(nums, k, 1, current, result);
    return result;
}

vector<string> node_ids(vector<Node*>& nodes) {
    vector<string> result;
    for (int i = 0; i < nodes.size(); i++) {
        result.push_back(nodes[i]->name);
    }
    return result;
}

void print_vector(vector<string>& values) {
    int i;
    cout << "[";
    for (i = 0; i < values.size() - 1; i++) {
        cout << values[i] << ", ";
    }
    cout << values[i] << "] ";
}

void process_data(vector<vector<Node*>> nodes_array, double& global_min, vector<string>& global_data) {
    double min_aic = INFINITY;
    vector<string> min_data;
    vector<Node*> combination;
    double aic; vector<string> array_id;
    for (int i = 0; i < nodes_array.size(); i++) {
        combination = nodes_array[i];
        aic = calculate_AIC(combination);
        array_id = node_ids(combination);
        // print_vector(array_id);
        // cout << setprecision(20) << "AIC value: " << aic << endl;
        if (aic < min_aic) {
            min_aic = aic;
            min_data = array_id;
        }
    }
    if (min_aic < global_min) {
        global_min = min_aic;
        global_data = min_data;
    }
    // cout << "" << endl;
    print_vector(min_data);
    cout << setprecision(20) << min_aic << endl;
}

void inorder(Node* root) {
    if (root == nullptr) {
        return;
    }
    inorder(root->left);
    cout << root->id << " ";
    inorder(root->right);
}

pair<double, Node*> forwardPartition(Node* currNode, vector<Node*> cutNodes, pair<double, Node*> bestPair) {
    if (currNode->left != nullptr) {
        bestPair = forwardPartition(currNode->left, cutNodes, bestPair);
    }
    if (currNode->right != nullptr) {
        bestPair = forwardPartition(currNode->right, cutNodes, bestPair);
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

pair<double, vector<string>> forwardRecursive(vector<Node*>& node_array) {
    const pair<double, Node*> defaultPair = {INFINITY, new Node(999, 999, "Default")};

    vector<Node*> cuts = {node_array[0]};
    vector<double> AICs;
    
    AICs.push_back(calculate_AIC(cuts));

    pair<double, Node*> tempPair;
    int numNodes = node_array.size();
    double tempAIC;
    pair<double, Node*> bestCutoff = {AICs.front(), cuts.front()};
    for (int i = 1; i < numNodes; i++) {
        tempPair = defaultPair; //becuase nothing has been done yet for this iteration
        tempPair = forwardPartition(node_array[0], cuts, tempPair);
        tempAIC = tempPair.first;
        AICs.push_back(tempAIC);
        cuts.push_back(tempPair.second);

        if (tempAIC < bestCutoff.first) {
            bestCutoff.first = tempAIC;
            bestCutoff.second = tempPair.second;
        }
    }

    //TODO: Will be commented out later
    int j = 0;
    for (Node* cutNode : cuts) {
        cout<<j<<" cuts ("<<j+1<<" rates): " << cutNode->name << ": " << setprecision(20) << AICs.at(j) << endl;
        j++;
    }

   vector<string> retVec;
    for (Node* cutNode : cuts) {
        retVec.push_back(cutNode->name);
        if (cutNode->name == bestCutoff.second->name) {
            break;
        }
    }
    pair<double, vector<string>> retPair = {bestCutoff.first, retVec};
    return retPair;
}


pair<double, Node*> backwardMerge(Node* currNode, Node* root, vector<Node*> cutNodes, pair<double, Node*> bestPair) {
     if (currNode->left != nullptr) {
        bestPair = backwardMerge(currNode->left, root, cutNodes, bestPair);
    }
    if (currNode->right != nullptr) {
        bestPair = backwardMerge(currNode->right, root, cutNodes, bestPair);
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


pair<double, vector<string>> backwardRecursive(vector<Node*>& node_array) {
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
        tempPair = backwardMerge(node_array[0], node_array[0], cuts, tempPair);
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

    //TODO: Will be commented out later
    int m = 0;
    for (Node* uncutNode: uncuts) {
        cout<<numCuts-m-1<<" cuts ("<<numCuts-m<<" rates): " << uncutNode->name << ": " << AICs.at(m) << endl;
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
    pair<double, vector<string>> retPair = {bestMerge.first, retVec};
    return retPair;
}

// [[Rcpp::export]]
void mainFunc(string newick_string) {

    vector<Node*> node_array;
    double global_min = INFINITY;
    vector<string> global_data;
    vector<vector<Node*>> result;
    // unordered_map<int, vector<vector<Node*>>> final_result;
    double time_elapsed;
    int size;
    
    clock_t start = clock();
    
    parse_newick(newick_string, node_array);
    node_array.erase(node_array.begin() + 1, node_array.begin() + 2);
    size = node_array.size();

    pair<double, vector<string>> forwardResult;
    cout << "Using Forward Algorithm:" << endl;
    cout << "" << endl;
    forwardResult = forwardRecursive(node_array);
    print_vector(forwardResult.second);
    cout << "is the " << forwardResult.second.size() << " partition array with the minimum AIC value: " << setprecision(20) << forwardResult.first << endl;
    cout << "" << endl;

    pair<double, vector<string>> backwardResult;
    cout << "Using Backward Algorithm:" << endl;
    cout << "" << endl;
    backwardResult = backwardRecursive(node_array);
    print_vector(backwardResult.second);
    cout << "is the " << backwardResult.second.size() << " partition array with the minimum AIC value: " << setprecision(20) << backwardResult.first << endl;
    cout << "" << endl;
    
    cout << "Using Bruteforce Algorithm:" << endl;
    cout << "" << endl;
    for (int i = 0; i < size - 1; i++) {
        result = combinations(node_array, i);
        // final_result[i+1] = result;
        // cout << "" << endl;
        // cout << "Total possible combinations for " << (i+1) << " partitions is: " << result.size() << endl;
        // cout << "" << endl;
        process_data(result, global_min, global_data);
    }
    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    print_vector(global_data);
    cout << "is the " << global_data.size() << " partition array with the minimum AIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed running 3 algorithms is: " << time_elapsed << endl;
    cout << "" << endl;
}