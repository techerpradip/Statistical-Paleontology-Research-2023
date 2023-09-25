#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <iomanip>
#include <chrono>
#include <thread>
#include "pthread.h"

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

struct ThreadData {
    vector<Node*>& tree_nodes;
    vector<vector<Node*>>& final_result;
    int partition_count;
    ThreadData(vector<Node*>& tree_nodes, vector<vector<Node*>>& final_result, int partition_count): tree_nodes(tree_nodes), final_result(final_result), partition_count(partition_count) {}
};

double global_min = INFINITY;
vector<Node*> global_data;
pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;


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
            split_string(current_string, name, rate, ':');
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

void generate_combinations(vector<Node*>& nums, vector<vector<Node*>>& final_result, int k, int start, int end, int index, vector<Node*>& current, double& min_aic) {
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
        generate_combinations(nums, final_result, k, i+1, end, index+1, current, min_aic);
    }
}


void* process_tree_parallel(void* arg) {
    ThreadData current_thread = *(ThreadData*) arg;
    int partition_count = current_thread.partition_count;
    vector<Node*>& nodes_array = current_thread.tree_nodes;
    vector<vector<Node*>>& final_result = current_thread.final_result;

    vector<Node*> current(partition_count, nullptr);
    double min_aic = INFINITY;
    current[0] = nodes_array[0];
    generate_combinations(nodes_array, final_result, partition_count, 1, nodes_array.size() - 1, 1, current, min_aic);

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

    delete defaultPair.second;

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

    delete defaultPair.second;
    delete uncuts[0];
    
    pair<double, vector<string>> retPair = {bestMerge.first, retVec};
    return retPair;
}

void delete_tree_nodes(Node* node) {
    if (node != nullptr) {
        delete_tree_nodes(node->left);
        delete_tree_nodes(node->right);
        delete node;
    }
}

void delete_thread_information(vector<ThreadData*>& thread_information) {
    for (int i = 0; i < thread_information.size(); i++) {
        delete thread_information[i];
    }
}

int main() {

    // string newick_string = "(Psarolepis_romeri:1,(Diabolepis_speratus:1.85,((Dipnorhynchus_kiandrensis:7.4,(Archaeonectes_pertusus:28.3,(Uranolophus_wyomingensis:1.6,(Speonesydrion_iani:0.8,(Jarvikia_arctica:36.5173913,(((Adololopas_moyasmithae:10.775,((Adelargo_schultzei:14.05,Chirodipterus_australis:3.25):6.1,(Chirodipterus_rhenanus:1.425,(Chirodipterus_wildungensis:3.25,Dipterus_cf_valenciennesi:3.25):4.675):1.425):1.425):10.31485507,(Barwickia_downunda:10.14492754,Dipterus_valenciennesi:4.444927536):4.444927536):4.444927536,(Pillararhynchus_longi:25.35217391,(((Gogodipterus_paddyensis:24.80434783,((Tarachomylax_oepiki:1.947826087,(Amadeodipterus_kencampbelli:0.9739130435,Stomiahykus_thlaodus:10.47391304):0.9739130435):0.9739130435,(Iowadipterus_halli:17.93913043,((Delatitia_breviceps:50.17391304,(Phaneropleuron_andersoni:34.69130435,((Orlovichthys_limnatis:34.32608696,(Howidipterus_donnae:16.84347826,(((Andreyevichthys_epitomus:16.88913043,Oervigia_nordica:16.88913043):16.88913043,(Grossipterus_crassus:22.79565217,(Fleurantia_denticulata:22.61304348,((Robinsondipterus_longi:16.275,(Asthenorhynchus_meemannae:10.85,(Holodipterus_elderae:5.425,Holodipterus_gogoensis:5.425):5.425):5.425):6.155434783,((Griphognathus_minutidens:14.46666667,(Griphognathus_sculpta:7.233333333,Griphognathus_whitei:7.233333333):7.233333333):7.78115942,(Rhynchodipterus_elginensis:32.86521739,(Jessenia_concentrica:0.1826086957,Soederberghia_groenlandica:21.8826087):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957,(Pentlandia_macroptera:8.330434783,Scaumenacia_curta:14.83043478):8.330434783):0.1826086957):0.1826086957):0.1826086957,(Holodipterus_santacrucensis:21.07439614,((Ganopristodus_splendens:55.8057971,(Megapleuron_zangerli:86.77149758,(Sagenodus_inaequalis:50.53719807,(((Eoctenodus_microsoma:2.634299517,Tranodis_castrensis:61.53429952):2.634299517,(Ctenodus_romeri:15.68429952,Straitonia_waterstoni:29.58429952):15.68429952):2.634299517,((Parasagenodus_sibiricus:11.33429952,(Gnathorhiza_serrata:20.55,((Beltanodus_ambilobensis:53.1125,(Namatozodia_pitikanta:45.525,(Ariguna_formosa:37.9375,(((Aphelodus_anapes:11.38125,Ceratodus_formosa:11.38125):11.38125,((Asiatoceratodus_sharovi:7.5875,Gosfordia_truncata:30.5875):7.5875,(Neoceratodus_forsteri:77.0875,(Mioceratodus_gregoryi:37,(Lepidosiren_paradoxa:22.4,Protopterus_annectens:9.5):9.5):86.5875):77.0875):7.5875):7.5875,(Archaeoceratodus_avus:26.675,Tellerodus_sturi:38.175):26.675):7.5875):7.5875):7.5875):12.3875,(Microceratodus_angolensis:63.9,(Palaeophichthys_parvulus:1.6,(Ptychoceratodus_serratus:45.525,(Paraceratodus_germaini:31.65,(Arganodus_atlantis:15.175,Ferganoceratodus_jurassicus:104.975):15.175):15.175):16.775):1.6):1.6):22.15):31.88429952):11.33429952,(Ceratodus_latissimus:76.93429952,Metaceratodus_wollastoni:159.4342995):76.93429952):11.33429952):2.634299517):2.634299517):2.634299517):2.634299517,(Nielsenia_nordica:14.62004831,Conchopoma_gadiforme:77.42004831):14.62004831):2.634299517):2.634299517):0.1826086957):0.1826086957):0.1826086957,(Rhinodipterus_secans:15.37826087,Rhinodipterus_ulrichi:8.87826087):8.87826087):0.1826086957):0.1826086957):0.1826086957):0.1826086957,(Palaeodaphus_insignis:12.49347826,Sunwapta_grandiceps:23.29347826):12.49347826):0.1826086957,(Melanognathus_canadensis:1.734782609,Sorbitorhynchus_deleaskitus:1.734782609):1.734782609):0.1826086957):0.1826086957):0.1826086957):0.9826086957):0.8):0.8):0.8):0.8,(Westollrhynchus_lehmanni:2,(Ichnomylax_kurnai:3.36,(Dipnorhynchus_sussmilchi:2.52,(Chirodipterus_onawwayensis:16.88,(Dipnorhynch_cathlesae:0.84,Dipnorhynchus_kurikae:0.84):0.84):0.84):0.84):2.84):2):2.65):1.85);";
    string newick_string = "((1:5,((2:10,(3:6,4:5)25:5)24:3,(5:8,((6:2,7:7)28:4,(8:1,9:9)29:5)27:7)26:8)23:6)22:2,((10:7,(11:10,(12:2,(13:10,14:8)34:7)33:7)32:8)31:8,(15:4,(16:3,(17:8,(18:3,(19:5,20:8)39:1)38:7)37:2)36:4)35:5)30:10)21;";
    // string newick_string = "((1:1,(2:3,(3:10,4:5)34:2)33:2)32:10,(5:2,((6:9,7:3)37:1,((8:9,((9:6,(10:7,(11:2,12:5)43:6)42:6)41:7,((13:3,14:10)45:2,(15:2,16:1)46:10)44:1)40:10)39:9,(17:4,((18:3,19:6)49:7,((20:2,(21:5,22:9)52:2)51:1,(23:6,((24:4,25:6)55:10,((26:7,27:1)57:5,(28:5,(29:6,30:8)59:8)58:10)56:4)54:2)53:7)50:1)48:5)47:3)38:1)36:7)35:4)31;";

    vector<Node*> node_array;
    int size;

    pthread_t* threads;
    
    const auto start = chrono::steady_clock::now();
    
    Node* root = parse_newick(newick_string, node_array);
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

    threads = (pthread_t*) malloc(sizeof(pthread_t) * size);

    vector<vector<Node*>> final_result(size, {node_array[0]});

    ThreadData current_thread(node_array, final_result, 0);
    vector<ThreadData> threads_information(size, current_thread);

    /*
     * Idea:
     *      Parallelise both combinations and pthread_create by creating a parallel wrapping function,
     *      maybe create a void* threads_wrapper(arg *)
    */

    // const int processor_count = thread::hardware_concurrency();

    for (int i = 1; i <= size; i++) {
        threads_information[i-1].partition_count = i;
        int create_result = pthread_create(&threads[i-1], nullptr, process_tree_parallel, &threads_information[i-1]);
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

    // for (int i = 2; i <= 3; i++) {
    //     threads_information[i-2].partition_count = i;
    //     int create_result = pthread_create(&threads[i-2], nullptr, process_tree_parallel, &threads_information[i-2]);
    //     if (create_result != 0) {
    //         cerr << "Error creating thread!" << endl;
    //         exit(-1);
    //     }
    // }

    // for (int i = 0; i < 2; i++) {
    //     int join_result = pthread_join(threads[i], nullptr);
    //     if (join_result != 0) {
    //         cerr << "Error joining thread: " << threads[i] << endl;
    //         exit(-1);
    //     }
    // }

    const auto end = chrono::steady_clock::now();

    cout << "" << endl;
    cout << "Among all the possible combinations in total," << endl;
    vector<string> global_min_ids = node_ids(global_data);
    print_vector(global_min_ids);
    cout << "is the " << global_data.size() << " partition array with the minimum AIC value: " << setprecision(20) << global_min << endl;
    cout << "" << endl;

    const chrono::duration<double> time_elapsed = end - start;
    cout << "Total time elapsed running 3 algorithms is: " << (time_elapsed.count())/60 << " minutes." << endl;
    cout << "" << endl;
    
    free(threads);
    threads = nullptr;

    delete_tree_nodes(root);

    // delete_thread_information(threads_information);

    pthread_mutex_destroy(&my_mutex);

    return 0;
}