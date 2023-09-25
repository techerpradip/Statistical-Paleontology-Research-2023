#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;

struct Node {
    int id;
    int rate;
    string name;
    Node* left = nullptr;
    Node* right = nullptr;
    Node(int id, int rate, string name): id(id), rate(rate), name(name) {}
};

void create_tree(vector<string>& arr_names, vector<int>& arr_rates, vector<Node*>& node_array) {
    int size = arr_names.size();
    vector<Node*> temp_vector;
    Node* temp = nullptr;
    int current_id = 1;
    for (int i = 0; i < size; i++) {
        if (arr_names[i] == "NULL") {
            temp_vector.push_back(nullptr);
        } else {
            temp = new Node(current_id, arr_rates[i], arr_names[i]);
            temp_vector.push_back(temp);
            node_array.push_back(temp);
            current_id++;
        }
    }
    
    for (int i = 0, j = 1; j < size; i++) {
        if (!temp_vector[i]) {
            continue;
        }
        temp_vector[i]->left = temp_vector[j++];
        if (j < size) {
            temp_vector[i]->right = temp_vector[j++];
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
    cout << "" << endl;
    print_vector(min_data);
    cout << "is the " << min_data.size() << " partition array with the minimum AIC value: " << setprecision(20) << min_aic << endl;
    cout << "" << endl;
    cout << "--------------------------------------------------" << endl;
}

void inorder(Node* root) {
    if (root == nullptr) {
        return;
    }
    inorder(root->left);
    cout << root->id << " ";
    inorder(root->right);
}


int main() {
    vector<string> arr_names = {"Jan", "Feb", "Mar", "NULL", "NULL", "Apr", "May", "Aug", "June", "Sep", "July"};
    vector<int> arr_rates = {0,7,3,0,0,15,5,4,10,3,2};

    // vector<string> arr_names;
    // vector<int> arr_rates;

    // for (int i = 0; i < 171; i++) {
    //     arr_names.push_back("Hello");
    //     arr_rates.push_back(10);
    // }

    vector<Node*> node_array;
    double global_min = INFINITY;
    vector<string> global_data;
    vector<vector<Node*>> result;
    Node* root = nullptr;
    // unordered_map<int, vector<vector<Node*>>> final_result;
    double time_elapsed;
    int size;
    
    clock_t start = clock();
    
    create_tree(arr_names, arr_rates, node_array);
    size = node_array.size();
    // root = parse_newick(newick_string, node_array);
    
    for (int i = 0; i < size; i++) {
        result = combinations(node_array, i);
        // final_result[i+1] = result;
        cout << "" << endl;
        cout << "Total possible combinations for " << (i+1) << " partitions is: " << result.size() << endl;
        cout << "" << endl;
        process_data(result, global_min, global_data);
    }
    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    print_vector(global_data);
    cout << "is the " << global_data.size() << " partition array with the minimum AIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed is: " << time_elapsed << endl;
    cout << "" << endl;

    return 0;
}