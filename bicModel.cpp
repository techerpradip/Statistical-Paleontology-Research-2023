#include <iostream>
#include <unordered_set>
#include <unordered_map> 
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;

struct Node {
    int id;
    int rate;
    int rcount = 0;
    int lcount = 0;
    Node* left = nullptr;
    Node* right = nullptr;

    Node(int id, int rate): id(id), rate(rate) {}
};

bool isPBT(int count) {
    count = count + 1;
    while (count % 2 == 0) {
        count = count / 2;
    }
    return (count == 1);
}

Node* insert(Node* root, int id, int rate, vector<Node*>& node_array) {
    if (root == nullptr) {
        Node* newNode = new Node(id, rate);
        node_array.push_back(newNode);
        return newNode;
    }

    if (root->rcount == root->lcount) {
        root->left = insert(root->left, id, rate, node_array);
        root->lcount += 1;
    } else if (root->rcount < root->lcount) {
        if (isPBT(root->lcount)) {
            root->right = insert(root->right, id, rate, node_array);
            root->rcount += 1;
        } else {
            root->left = insert(root->left, id, rate, node_array);
            root->lcount += 1;
        }
    }
    
    return root;
}

int sum_of_rates(Node* node, unordered_set<Node*>& partition_set, int i) {
    if (node == nullptr) {
        return 0;
    }
    if (partition_set.find(node) != partition_set.end() && i == 1) {
        return 0;
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

double calculate_BIC(vector<Node*>& combination) {
    Node* node;
    int sum; int edges; double average; double bic_value;
    double likelihood = 0;
    double correction;
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
    correction = (2 * k * (k + 1))/(12 - k - 1);
    bic_value = 2 * (k - likelihood) + correction;
    // bic_value =  (log(8) * k) - (2 * likelihood); // COME BACK LATER
    return bic_value;
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

vector<int> node_ids(vector<Node*>& nodes) {
    vector<int> result;
    for (int i = 0; i < nodes.size(); i++) {
        result.push_back(nodes[i]->id);
    }
    return result;
}

void print_vector(vector<int>& values) {
    int i;
    cout << "[";
    for (i = 0; i < values.size() - 1; i++) {
        cout << values[i] << ", ";
    }
    cout << values[i] << "] ";
}

void process_data(vector<vector<Node*>> nodes_array, double& global_min, vector<int>& global_data) {
    double min_bic = INFINITY;
    vector<int> min_data;
    vector<Node*> combination;
    double bic; vector<int> array_id;
    for (int i = 0; i < nodes_array.size(); i++) {
        combination = nodes_array[i];
        bic = calculate_BIC(combination);
        array_id = node_ids(combination);
        // print_vector(array_id);
        // cout << setprecision(20) << "BIC value: " << bic << endl;
        if (bic < min_bic) {
            min_bic = bic;
            min_data = array_id;
        }
    }
    if (min_bic < global_min) {
        global_min = min_bic;
        global_data = min_data;
    }
    cout << "" << endl;
    print_vector(min_data);
    cout << "is the " << min_data.size() << " partition array with the minimum BIC value: " << setprecision(20) << min_bic << endl;
    cout << "" << endl;
    cout << "--------------------------------------------------" << endl;
}

// void delete_nodes(Node* node) {
//     if (node != nullptr) {
//         delete_nodes(node->left);
//         delete_nodes(node->right);
//         delete node;
//     }
// }

int main() {
    // vector<int> arr_ids = {1, 3, 2, 4, 7, 6, 5, 9, 8};
    // vector<int> arr_rates = {0, 5, 2, 6, 1, 3, 0, 7, 8};

    vector<int> arr_ids;
    vector<int> arr_rates;
    for (int i = 0; i < 13; i++) {
        arr_ids.push_back(i+1);
        arr_rates.push_back(i+1);
        // arr_rates.push_back(rand() % 10);
    }

    // vector<int> arr_ids = {1, 3, 2};
    // vector<int> arr_rates = {0, 5, 2};
    int size = arr_ids.size();
    Node* root = nullptr;
    vector<Node*> node_array;
    double global_min = INFINITY;
    vector<int> global_data;
    vector<vector<Node*>> result;
    unordered_map<int, vector<vector<Node*>>> final_result;
    double time_elapsed;
    
    clock_t start = clock();
    for (int i = 0; i < size; i++) {
        root = insert(root, arr_ids[i], arr_rates[i], node_array);
    }
    
    for (int i = 0; i < size; i++) {
        result = combinations(node_array, i);
        final_result[i+1] = result;
        cout << "" << endl;
        cout << "Total possible combinations for " << (i+1) << " partitions is: " << result.size() << endl;
        cout << "" << endl;
        process_data(result, global_min, global_data);
    }
    clock_t end = clock();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    print_vector(global_data);
    cout << "is the " << global_data.size() << " partition array with the minimum BIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    time_elapsed = double(end - start) / double(CLOCKS_PER_SEC);
    cout << setprecision(10) << "Total time elapsed is: " << time_elapsed << endl;
    cout << "" << endl;

    return 0;
}