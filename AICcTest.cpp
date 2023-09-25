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
        name = tree_name + " ";
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
    Node* root = nullptr;
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
            new_node = new Node(current_id, rate, name + to_string(current_id));
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
                new_node = new Node(current_id, rate, name + to_string(current_id));
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
                new_node = new Node(current_id, rate, name + to_string(current_id));
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
    correction = (2 * k * (k + 1))/(8 - k - 1);
    aic_value = 2 * (k - likelihood) + correction;
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
        print_vector(array_id);
        cout << setprecision(20) << "AIC value: " << aic << endl;
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

    // string newick_string = "(Psarolepis_romeri:1,(Diabolepis_speratus:1.85,((Dipnorhynchus_kiandrensis:7.4,(Archaeonectes_pertusus:28.3,(Uranolophus_wyomingensis:1.6,(Speonesydrion_iani:0.8,(Jarvikia_arctica:36.5173913,(((Adololopas_moyasmithae:10.775,((Adelargo_schultzei:14.05,Chirodipterus_australis:3.25):6.1,(Chirodipterus_rhenanus:1.425,(Chirodipterus_wildungensis:3.25,Dipterus_cf_valenciennesi:3.25):4.675):1.425):1.425):10.31485507,(Barwickia_downunda:10.14492754,Dipterus_valenciennesi:4.444927536):4.444927536):4.444927536,(Pillararhynchus_longi:25.35217391,(((Gogodipterus_paddyensis:24.80434783,((Tarachomylax_oepiki:1.947826087,(Amadeodipterus_kencampbelli:0.9739130435,Stomiahykus_thlaodus:10.47391304):0.9739130435):0.9739130435,(Iowadipterus_halli:17.93913043,((Delatitia_breviceps:50.17391304,(Phaneropleuron_andersoni:34.69130435,((Orlovichthys_limnatis:34.32608696,(Howidipterus_donnae:16.84347826,(((Andreyevichthys_epitomus:16.88913043,Oervigia_nordica:16.88913043):16.88913043,(Grossipterus_crassus:22.79565217,(Fleurantia_denticulata:22.61304348,((Robinsondipterus_longi:16.275,(Asthenorhynchus_meemannae:10.85,(Holodipterus_elderae:5.425,Holodipterus_gogoensis:5.425):5.425):5.425):6.155434783,((Griphognathus_minutidens:14.46666667,(Griphognathus_sculpta:7.233333333,Griphognathus_whitei:7.233333333):7.233333333):7.78115942,(Rhynchodipterus_elginensis:32.86521739,(Jessenia_concentrica:0.1826086957,Soederberghia_groenlandica:21.8826087):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957,(Pentlandia_macroptera:8.330434783,Scaumenacia_curta:14.83043478):8.330434783):0.1826086957):0.1826086957):0.1826086957,(Holodipterus_santacrucensis:21.07439614,((Ganopristodus_splendens:55.8057971,(Megapleuron_zangerli:86.77149758,(Sagenodus_inaequalis:50.53719807,(((Eoctenodus_microsoma:2.634299517,Tranodis_castrensis:61.53429952):2.634299517,(Ctenodus_romeri:15.68429952,Straitonia_waterstoni:29.58429952):15.68429952):2.634299517,((Parasagenodus_sibiricus:11.33429952,(Gnathorhiza_serrata:20.55,((Beltanodus_ambilobensis:53.1125,(Namatozodia_pitikanta:45.525,(Ariguna_formosa:37.9375,(((Aphelodus_anapes:11.38125,Ceratodus_formosa:11.38125):11.38125,((Asiatoceratodus_sharovi:7.5875,Gosfordia_truncata:30.5875):7.5875,(Neoceratodus_forsteri:77.0875,(Mioceratodus_gregoryi:37,(Lepidosiren_paradoxa:22.4,Protopterus_annectens:9.5):9.5):86.5875):77.0875):7.5875):7.5875,(Archaeoceratodus_avus:26.675,Tellerodus_sturi:38.175):26.675):7.5875):7.5875):7.5875):12.3875,(Microceratodus_angolensis:63.9,(Palaeophichthys_parvulus:1.6,(Ptychoceratodus_serratus:45.525,(Paraceratodus_germaini:31.65,(Arganodus_atlantis:15.175,Ferganoceratodus_jurassicus:104.975):15.175):15.175):16.775):1.6):1.6):22.15):31.88429952):11.33429952,(Ceratodus_latissimus:76.93429952,Metaceratodus_wollastoni:159.4342995):76.93429952):11.33429952):2.634299517):2.634299517):2.634299517):2.634299517,(Nielsenia_nordica:14.62004831,Conchopoma_gadiforme:77.42004831):14.62004831):2.634299517):2.634299517):0.1826086957):0.1826086957):0.1826086957,(Rhinodipterus_secans:15.37826087,Rhinodipterus_ulrichi:8.87826087):8.87826087):0.1826086957):0.1826086957):0.1826086957):0.1826086957,(Palaeodaphus_insignis:12.49347826,Sunwapta_grandiceps:23.29347826):12.49347826):0.1826086957,(Melanognathus_canadensis:1.734782609,Sorbitorhynchus_deleaskitus:1.734782609):1.734782609):0.1826086957):0.1826086957):0.1826086957):0.9826086957):0.8):0.8):0.8):0.8,(Westollrhynchus_lehmanni:2,(Ichnomylax_kurnai:3.36,(Dipnorhynchus_sussmilchi:2.52,(Chirodipterus_onawwayensis:16.88,(Dipnorhynch_cathlesae:0.84,Dipnorhynchus_kurikae:0.84):0.84):0.84):0.84):2.84):2):2.65):1.85);";
    string newick_string = "(((July:2,Sep:3)May:5,(June:10,Aug:4)Apr:15)Mar:3,Feb:7)Jan;";

    vector<Node*> node_array;
    double global_min = INFINITY;
    vector<string> global_data;
    vector<vector<Node*>> result;
    Node* root = nullptr;
    // unordered_map<int, vector<vector<Node*>>> final_result;
    double time_elapsed;
    int size;
    
    clock_t start = clock();
    
    root = parse_newick(newick_string, node_array);
    size = node_array.size();
    
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