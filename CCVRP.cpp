

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <string>
#include <string.h>
#include <algorithm>
#include <random>
#include <iomanip>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <chrono>
#include <functional>
#include <unordered_map>
#include <cstdlib>
#include <numeric>


using namespace std;

#define FUNC_TIMER()  FunctionTimer __timer(__func__)

clock_t startTime, endTime;
int timeMax = 10;
int t = 0;

//string inFile = "CMT1.vrp";
//string instName = "CMT";
//string instName = "Li_";
//string instName = "F_";
string instName = "B";
string instanceName = "B-n68-k9";
//string instName = "Golden_";
int instNum = 1;
string instNumber = to_string(instNum);
string instType = ".vrp";
//string inFile = instName + "/" + instName + instNumber + instType;
string inFile = instName + "/" + instanceName + instType;

//string inst = instName + instNumber;
string inst = instanceName;
string vehicle_num_file = "Instance Vehicle Num/" + instName + ".txt";

int node_num;                                     // vertice number
int customer_num;                                 // customer number
vector<double> xCor;                                                   // x
vector<double> yCor;                                                   // y
vector<double> demand;                                                 // demand
vector<vector<double>> travelTime;                                     
vector<double> theta;                                                  
int load_capacity;
int vehicle_num;                                                    
double M = 10000000;                                               
int max_iteration = 100000;
int** Tabu;                                                        
int tabu_step = 1;                                                  
int s_r;                                                           

int op_num = 15;
vector<int> op_used_num(op_num, 0);                               
vector<int> op_effec_num(op_num, 0);                              

vector<vector<int>> sorted_neighbors;                             

struct routeInformation {
	vector<vector<int>> route;                                      
	vector<vector<double>> arrival_time;                            
	vector<int> used_capacity;                                      
	vector<double> single_cumulative_time;                          
};
routeInformation S;

struct evaluationFunction {
	double cumulative_time;                                         
	double extra_capacity;                                          
	double extra_vehicle;                                           
	double evaluation_obj;                                          
};
evaluationFunction F;

struct globalBest {
	vector<vector<int>> route;                                      
	double obj;                                                    
};
globalBest optimal;


struct penaltyParameter {
	double alpha_extra_capacity = 40;                               
	double alpha_extra_vehicle = 800;                              
	double extra_capacity_delta;                                    
	double extra_vehicle_delta;                                     
	double max_extra_capacity;
	double min_extra_capacity;
	double max_extra_vehicle;
	double min_extra_vehicle;
};
penaltyParameter penalty;

struct tempInformation {
	vector<int> single_route;                                       
	vector<double> single_arrival_time;                             
	double cumulative_time;
};

struct operatorPara {
	int last_node;              
	int next_node;              
	bool isInsert;              
	bool isChange;              
	double delta;              
	double tem_delta;           
	int select_last_node_index; 
	int select_next_node_index; 
	int select_route_index;     
	int select_node_index;      
	int node1_index;            
	int node2_index;            
	int node1;                  
	int node2;                  
	int route1_index;           
	int route2_index;           
};

struct parameters {
	double shake_s = 4;                                                 
	int max_shake_no_improve = 100;                                     
	double max_length_RCL_percent = 0.20;                                   
	int max_initial_LS = 1000;                                         
	int max_no_improve_new_S = 2000;                                    
	int omega = 3;                                                    
	int min_omega = 2;
	int max_omega = 6;                                               
	int crossover_no_improve_cycle = 5;                                     
} p;

// ====================== increment evaluation ======================




pair<double, double> calc_remove_delta(const routeInformation& S, int r_idx, int n_idx) {
	if (r_idx >= S.route.size() || n_idx <= 0 || n_idx >= S.route[r_idx].size() - 1) return { 0.0, 0.0 };

	int node = S.route[r_idx][n_idx];
	int prev = S.route[r_idx][n_idx - 1];
	int next = S.route[r_idx][n_idx + 1];

	double arrival_time_node = S.arrival_time[r_idx][n_idx];

	// time delta
	double delta_time = travelTime[prev][next] - (travelTime[prev][node] + travelTime[node][next]);

	
	int nodes_after = S.route[r_idx].size() - 1 - n_idx - 1; 
	

	int customers_after = 0;
	
	
	if (n_idx + 1 <= S.route[r_idx].size() - 2) {
		customers_after = (S.route[r_idx].size() - 2) - (n_idx + 1) + 1;
	}

	
	double delta_cum_time = -arrival_time_node + delta_time * customers_after;

	return { -(double)demand[node], delta_cum_time };
}


pair<double, double> calc_insert_delta(const routeInformation& S, int r_idx, int pos_idx, int node) {
	

	int prev = S.route[r_idx][pos_idx - 1];
	int next = S.route[r_idx][pos_idx]; 

	double arrival_time_prev = S.arrival_time[r_idx][pos_idx - 1];
	double arrival_time_node = arrival_time_prev + travelTime[prev][node];

	double delta_time = (travelTime[prev][node] + travelTime[node][next]) - travelTime[prev][next];

	int customers_after = 0;
	
	if (pos_idx <= S.route[r_idx].size() - 2) {
		customers_after = (S.route[r_idx].size() - 2) - pos_idx + 1;
	}

	
	double delta_cum_time = arrival_time_node + delta_time * customers_after;

	return { (double)demand[node], delta_cum_time };
}


double calc_penalty_delta(double current_load, double delta_load, int load_cap, double alpha) {
	double old_over = max(0.0, current_load - load_cap);
	double new_over = max(0.0, current_load + delta_load - load_cap);
	return (new_over - old_over) * alpha;
}

// ====================== key class ======================
class FunctionTimer {
public:
	explicit FunctionTimer(const char* func_name)
		: func_name_(func_name),
		start_(std::chrono::steady_clock::now()) {
	}

	~FunctionTimer() {
		auto end = std::chrono::steady_clock::now();
		auto us = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();

		getTotals()[func_name_] += us;      
		getCalls()[func_name_]++;           
	}

	
	static void PrintReport() {
		std::cout << "\n=== report ===\n";
		std::cout << std::left << std::setw(40) << "function name"
			<< std::setw(15) << "time (ms)"
			<< std::setw(12) << "number "
			<< std::setw(15) << "average (μs)" << "\n";
		std::cout << std::string(82, '-') << "\n";

		for (const auto& p : getTotals()) {           
			const std::string& name = p.first;
			long long total_us = p.second;

			double total_ms = total_us / 1000.0;
			int calls = getCalls()[name];
			double avg_us = calls ? static_cast<double>(total_us) / calls : 0.0;

			std::cout << std::left << std::setw(40) << name
				<< std::setw(15) << std::fixed << std::setprecision(2) << total_ms
				<< std::setw(12) << calls
				<< std::setw(15) << std::fixed << std::setprecision(1) << avg_us << "\n";
		}
	}

private:
	const char* func_name_;
	std::chrono::time_point<std::chrono::steady_clock> start_;

	static std::unordered_map<std::string, long long>& getTotals() {
		static std::unordered_map<std::string, long long> totals;
		return totals;
	}
	static std::unordered_map<std::string, int>& getCalls() {
		static std::unordered_map<std::string, int> calls;
		return calls;
	}
};


void readInitial_CMT(istream&);
void cumulativeTimeForRoute(const routeInformation&, evaluationFunction&);
void cumulativeTimeForSingleRoute(routeInformation&, const int);
void cumulativeTimeForSingleRoute(const vector<double>&, const vector<int>&, double&);
void extraPowerForRoute(const routeInformation&, evaluationFunction&);
void arrivalTimeForSingleRoute(vector<double>&, const vector<int>&);
void verify(const routeInformation&, const evaluationFunction&);
void verify(const routeInformation&);
void deltaEvaluation(tempInformation&, evaluationFunction&, const int);
void deltaEvaluation(tempInformation&, evaluationFunction&, const int, tempInformation&);
void deltaEvaluationForTwoSingle(tempInformation&, tempInformation&, evaluationFunction&, const int, const int);
void cumulativeTimeForSingleRoute(routeInformation&, int);
void printVerifyInformation(const routeInformation&, const operatorPara&, const evaluationFunction&);
void printVerifyInformationForTwoRoute(const routeInformation&, const operatorPara&, const evaluationFunction&);
void printVerifySingleRoute(const tempInformation&);
void calInformation(routeInformation&, evaluationFunction&);
void dEAX(const routeInformation&, const routeInformation&, routeInformation&, evaluationFunction&);
void update_route_arrival_time(routeInformation&, int);
// ================== operator ==================
void swap_star(routeInformation& S);
void node_arc_exchange(routeInformation& S);
void ejection_chain_mutation(routeInformation& S);
void route_reversal(routeInformation& S);
void init_neighbors();


double randomDouble() {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<double> dis(0.0, 1.0);

	return dis(gen);
}

int randomInt(int n) {
	if (n <= 0) return 0;
	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, n - 1);
	return dis(gen);
}

int randomInt(int lb, int ub) {
	if (lb > ub) {
		std::cerr << "Warning: lb > ub in randomInt(" << lb << ", " << ub << ")" << std::endl;
		return lb;
	}
	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(lb, ub);
	return dis(gen);
}

void readInitial_CMT(istream& ios) {
	double wait_time;
	ios >> node_num;
	customer_num = node_num - 1;
	ios >> load_capacity;
	ios >> wait_time;
	double a;

	for (int i = 0; i < node_num; i++) {

		ios >> a;
		ios >> a;
		xCor.push_back(a);
		ios >> a;
		yCor.push_back(a);
	}

	for (auto i = 0; i < node_num; i++) {
		ios >> a;
		ios >> a;
		demand.push_back(a);
	}

	
	if (Tabu != nullptr) {
		for (int i = 0; i < node_num; i++) {
			delete[] Tabu[i];
		}
		delete[] Tabu;
		Tabu = nullptr;
	}

	Tabu = new int* [node_num];
	for (auto i = 0; i < node_num; i++) {
		Tabu[i] = new int[node_num];
		for (auto j = 0; j < node_num; j++) {
			Tabu[i][j] = -1;
		}
	}

	if (node_num <= 50) {
		p.max_no_improve_new_S = 1000;
	}
	else if (node_num > 50 && node_num < 100) {
		p.max_no_improve_new_S = 5000;
	}
	else {
		p.max_no_improve_new_S = 10000;
	}

}

void readInitial_CMT(string) {
	ifstream FIC;
	FIC.open(inFile);
	if (FIC.fail()) {
		std::cout << "### Erreur open, File_Name: " << inFile << endl;
		getchar();
		exit(0);
	}
	if (FIC.eof()) {
		std::cout << "### Error open, File_Name: " << inFile << endl;
		getchar();
		exit(0);
	}
	readInitial_CMT(FIC);

}

void readVRP(const string& filename) {
	ifstream fin(filename);
	if (!fin.is_open()) {
		cerr << "Cannot open file: " << filename << endl;
		exit(1);
	}

	string line;
	string section = "";

	int dimension = 0;

	
	xCor.clear();
	yCor.clear();
	demand.clear();

	while (getline(fin, line)) {

		if (line.find("DIMENSION") != string::npos) {
			sscanf(line.c_str(), "DIMENSION : %d", &dimension);
			node_num = dimension;
			customer_num = node_num - 1;
		}

		else if (line.find("CAPACITY") != string::npos) {
			sscanf(line.c_str(), "CAPACITY : %d", &load_capacity);
		}

		else if (line.find("NODE_COORD_SECTION") != string::npos) {
			section = "NODE";
			xCor.resize(node_num);
			yCor.resize(node_num);
		}

		else if (line.find("DEMAND_SECTION") != string::npos) {
			section = "DEMAND";
			demand.resize(node_num);
		}

		else if (line.find("DEPOT_SECTION") != string::npos) {
			section = "DEPOT";
		}

		else if (line.find("EOF") != string::npos) {
			break;
		}

		else {
			if (section == "NODE") {
				int id;
				double x, y;
				if (sscanf(line.c_str(), "%d %lf %lf", &id, &x, &y) == 3) {
					xCor[id - 1] = x;
					yCor[id - 1] = y;
				}
			}

			else if (section == "DEMAND") {
				int id;
				double d;
				if (sscanf(line.c_str(), "%d %lf", &id, &d) == 2) {
					demand[id - 1] = d;
				}
			}

			else if (section == "DEPOT") {

			}
		}
	}

	fin.close();

	// ================= initial Tabu =================
	if (Tabu != nullptr) {
		for (int i = 0; i < node_num; i++) {
			delete[] Tabu[i];
		}
		delete[] Tabu;
	}

	Tabu = new int* [node_num];
	for (int i = 0; i < node_num; i++) {
		Tabu[i] = new int[node_num];
		for (int j = 0; j < node_num; j++) {
			Tabu[i][j] = -1;
		}
	}

	// ================= adjust parameter =================
	if (node_num <= 50) {
		p.max_no_improve_new_S = 1000;
	}
	else if (node_num < 100) {
		p.max_no_improve_new_S = 5000;
	}
	else {
		p.max_no_improve_new_S = 10000;
	}

	

	cout << "VRP file loaded: " << filename << endl;
	cout << "Nodes: " << node_num << ", Capacity: " << load_capacity << endl;
}

void getVehicleNum(const string& filename, const string& target) {
	ifstream fin(filename);
	if (!fin.is_open()) {
		cerr << "Cannot open file!" << endl;
		exit(1);
	}

	string name;
	int num;
	bool isFind = 0;

	while (fin >> name >> num) {
		if (name == target) {
			fin.close();
			vehicle_num = num;
			isFind = 1;
		}
	}

	fin.close();

	
	if (!isFind) {
		cerr << "Instance not found: " << target << endl;
		exit(1);
	}

}

void calDistance() {
	travelTime.resize(node_num);
	double rawTime;
	for (int i = 0; i < node_num; i++) {
		travelTime[i].reserve(node_num);
		for (int j = 0; j < node_num; j++) {
			rawTime = hypot(xCor[i] - xCor[j], yCor[i] - yCor[j]);
			//travelTime[i].push_back(round(rawTime * 10000) / 10000.0);
			//travelTime[i].push_back(round(rawTime * 100000) / 100000.0);
			travelTime[i].push_back(round(rawTime));
		}
	}
}

void GRASP_initialization(routeInformation& S) {
	vector<int> U;
	for (int i = 0; i < customer_num; i++) {
		U.push_back(i + 1);
	}

	for (auto i = 0; i < vehicle_num; i++) {
		S.route.emplace_back();
		S.route.back().push_back(0);
		S.route.back().push_back(0);


		S.used_capacity.push_back(0);
		S.single_cumulative_time.push_back(0);

	}

	vector<int> RCL;
	int max_length_RCL = p.max_length_RCL_percent * customer_num;

	vector<double> delta;                  
	vector<int> location_route;            
	//vector<int> tem_used_capacity;         
	int count;
	int node;
	int last_node;
	int last_index;
	int next_index;
	int next_node;
	int RCL_length;
	double tem_delta;
	int select_node;
	int select_node_index;
	int select_route;

	while (U.size()) {
		RCL.clear();
		delta.clear();
		location_route.clear();
		count = 0;

		for (auto i = 0; i < U.size(); i++) {
			node = U[i];
			for (auto j = 0; j < vehicle_num; j++) {
				last_node = S.route[j].back();                       
				last_index = S.route[j].size() - 1;                  
				next_index = S.route[j].size() - 2;                  
				next_node = S.route[j][next_index];                  
				tem_delta = travelTime[last_node][node] + travelTime[node][next_node] - travelTime[last_node][next_node];
				if (RCL.size() <= count) {
					RCL.push_back(node);
					delta.push_back(tem_delta);
					location_route.push_back(j);
					//tem_used_capacity.push_back(S.used_capacity[j] + demand[node]);
					count++;
				}
				else {
					if (tem_delta < delta[j]) {
						delta[j] = tem_delta;
						RCL[j] = node;
						location_route[j] = j;
						//tem_used_capacity[j] = S.used_capacity[j] + demand[node];
						count++;
					}
				}

				if (count >= vehicle_num) {
					count = 0;
				}
			}
		}

		
		select_node_index = randomInt(RCL.size());
		select_node = RCL[select_node_index];
		select_route = location_route[select_node_index];
		next_index = S.route[select_route].size() - 1;                        
		S.route[select_route].insert(S.route[select_route].begin() + next_index, select_node);
		S.used_capacity[select_route] += demand[select_node];

		for (auto k = 0; k < U.size(); k++) {
			if (U[k] == select_node) {
				U.erase(U.begin() + k);
				break;
			}
		}
	}

	for (auto i = 0; i < S.route.size(); i++) {
		S.arrival_time.emplace_back();                                 
		S.arrival_time.back().reserve(node_num);                      
	}

	int node1;
	int node2;
	double temp_arrival_time;
	for (auto i = 0; i < S.route.size(); i++) {
		temp_arrival_time = 0;
		S.arrival_time[i].push_back(0);
		for (auto j = 0; j < S.route[i].size() - 1; j++) {
			node1 = S.route[i][j];
			node2 = S.route[i][j + 1];
			temp_arrival_time += travelTime[node1][node2];
			S.arrival_time[i].push_back(temp_arrival_time);
		}
	}



	
	cumulativeTimeForRoute(S, F);

	for (auto i = 0; i < S.route.size(); i++) {
		cumulativeTimeForSingleRoute(S, i);
	}

	F.extra_capacity = 0;
	for (auto i = 0; i < S.route.size(); i++) {
		F.extra_capacity += max(0, S.used_capacity[i] - load_capacity);
	}

	F.extra_vehicle = max(0, static_cast<int>(S.route.size()) - vehicle_num);

	optimal.route = S.route;
	optimal.obj = F.cumulative_time;

	F.evaluation_obj = F.cumulative_time + penalty.alpha_extra_capacity * F.extra_capacity + penalty.alpha_extra_vehicle * F.extra_vehicle;

}


void calInformation(routeInformation& S, evaluationFunction& F) {
	S.used_capacity.clear();
	int tem_capacity;
	for (auto i = 0; i < S.route.size(); i++) {
		tem_capacity = 0;
		for (auto j = 1; j < S.route[i].size() - 1; j++) {
			tem_capacity += demand[S.route[i][j]];
		}
		S.used_capacity.push_back(tem_capacity);
	}

	S.arrival_time.clear();
	for (auto i = 0; i < S.route.size(); i++) {
		S.arrival_time.emplace_back();                                 
		S.arrival_time.back().reserve(node_num);                       
	}

	int node1;
	int node2;
	double temp_arrival_time;
	for (auto i = 0; i < S.route.size(); i++) {
		temp_arrival_time = 0;
		S.arrival_time[i].push_back(0);
		for (auto j = 0; j < S.route[i].size() - 1; j++) {
			node1 = S.route[i][j];
			node2 = S.route[i][j + 1];
			temp_arrival_time += travelTime[node1][node2];
			S.arrival_time[i].push_back(temp_arrival_time);
		}
	}

	
	cumulativeTimeForRoute(S, F);

	S.single_cumulative_time.clear();
	for (auto i = 0; i < S.route.size(); i++) {
		S.single_cumulative_time.push_back(0);
	}

	for (auto i = 0; i < S.route.size(); i++) {
		cumulativeTimeForSingleRoute(S, i);
	}

	F.extra_capacity = 0;
	for (auto i = 0; i < S.route.size(); i++) {
		F.extra_capacity += max(0, S.used_capacity[i] - load_capacity);
	}

	F.extra_vehicle = max(0, static_cast<int>(S.route.size()) - vehicle_num);

	F.evaluation_obj = F.cumulative_time + penalty.alpha_extra_capacity * F.extra_capacity + penalty.alpha_extra_vehicle * F.extra_vehicle;

	/*cout << "The information of the route: " << endl;
	for (auto i = 0; i < S.route.size(); i++) {
		for (auto j = 0; j < S.route[i].size(); j++) {
			cout << S.route[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	cout << "The information of the best route: " << endl;
	for (auto i = 0; i < optimal.route.size(); i++) {
		for (auto j = 0; j < optimal.route[i].size(); j++) {
			cout << optimal.route[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	cout << "The information of the arrival time: " << endl;
	for (auto i = 0; i < S.arrival_time.size(); i++) {
		for (auto j = 0; j < S.arrival_time[i].size(); j++) {
			cout << S.arrival_time[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	cout << "The information of whether the BSS is used: " << endl;
	for (auto i = 0; i < S.isBSS.size(); i++) {
		for (auto j = 0; j < S.isBSS[i].size(); j++) {
			cout << S.isBSS[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	cout << "The cumulative waiting time of the whole route: " << F.cumulative_time << endl;
	cout << "The best cumulative waiting time of the whole route: " << optimal.obj << endl;
	cout << endl;
	cout << "The cumulative waiting time of each route is: ";
	for (auto i = 0; i < S.route.size(); i++) {
		cout << S.single_cumulative_time[i] << " ";
	}
	cout << endl;
	cout << endl;

	cout << "The cumulative extra power of the whole route: " << F.extra_power << endl;
	cout << endl;
	cout << "The extra power of each route is: ";
	for (auto i = 0; i < S.route.size(); i++) {
		cout << S.single_extra_power[i] << " ";
	}
	cout << endl;

	cout << "The number of extra vehicle is: " << F.extra_vehicle << endl;
	cout << endl;

	cout << "The evaluation function is " << F.evaluation_obj << endl;

	double distance = 0;
	double single_distance = 0;
	for (auto i = 0; i < S.route.size(); i++) {
		single_distance = 0;
		for (auto j = 0; j < S.route[i].size() - 1; j++) {
			node1 = S.route[i][j];
			node2 = S.route[i][j + 1];
			distance += travelTime[node1][node2];
			single_distance += travelTime[node1][node2];
		}
		cout << "Route " << i << " 's distance is " << single_distance << endl;
	}
	cout << "The distance of the whole route is " << distance << endl;*/
}


void cumulativeTimeForRoute(const routeInformation& S, evaluationFunction& F) {
	double temp_cumulative_time = 0;
	auto len_S_route = S.route.size();
	for (auto i = 0; i < len_S_route; i++) {
		auto len_S_route_i = S.route[i].size() - 1;
		for (auto j = 1; j < len_S_route_i; j++) {
			temp_cumulative_time += S.arrival_time[i][j];
		}
	}
	F.cumulative_time = temp_cumulative_time;
}


void cumulativeTimeForSingleRoute(routeInformation& S, const int i) {
	double temp_single_cumulative_time = 0;
	auto len_S_route_i = S.route[i].size() - 1;
	for (auto j = 1; j < len_S_route_i; j++) {
		temp_single_cumulative_time += S.arrival_time[i][j];
	}
	S.single_cumulative_time[i] = temp_single_cumulative_time;
	
}


void cumulativeTimeForSingleRoute(const vector<double>& arrival_time, const vector<int>& route, double& single_cumulative_time) {
	double temp_single_cumulative_time = 0;
	auto len_route = route.size() - 1;
	for (auto j = 1; j < len_route; j++) {
		temp_single_cumulative_time += arrival_time[j];
	}
	single_cumulative_time = temp_single_cumulative_time;
	//cout << single_cumulative_time << endl;
}



void arrivalTimeForSingleRoute(vector<double>& arrival_time, const vector<int>& single_route) {
	int node1;
	int node2;
	double temp_arrival_time = 0;
	arrival_time.clear();
	arrival_time.push_back(0);
	for (auto j = 0; j < single_route.size() - 1; j++) {
		node1 = single_route[j];
		node2 = single_route[j + 1];
		temp_arrival_time += travelTime[node1][node2];
		arrival_time.push_back(temp_arrival_time);
	}

}


void verify(const routeInformation& S) {
	
	int tem_capacity = 0;
	bool isCapacity = true;                           
	for (auto i = 0; i < S.route.size(); i++) {
		tem_capacity = 0;
		for (auto j = 0; j < S.route[i].size(); j++) {
			tem_capacity += demand[S.route[i][j]];
		}
		if (tem_capacity > load_capacity) {
			isCapacity = false;
			std::cout << "The CAPACITY constraint for route " << i << " has been violated with " << (tem_capacity - load_capacity) << endl;
		}
	}


	if (!isCapacity) {
		if (!isCapacity) {
			std::cout << "△ Warning! The CAPACITY constraint is violated and the details are above ↑" << endl;
		}
	}

}


void verify(const routeInformation& S, const evaluationFunction& F) {
	
	bool isCapacity = true;
	double tem_capacity = 0;
	for (auto i = 0; i < S.route.size(); i++) {
		tem_capacity = 0;
		for (auto j = 0; j < S.route[i].size(); j++) {
			tem_capacity += demand[S.route[i][j]];
		}
		if (abs(tem_capacity - S.used_capacity[i]) > 0.01) {
			isCapacity = false;
			std::cout << "△ Warning! The CAPACITY for route " << i << " is wrongly calculated! Original: " << S.used_capacity[i] << " Validation: " << tem_capacity << endl;
		}
	}

	
	bool isExtraCapacity = true;
	double tem_extra_capacity = 0;
	for (auto i = 0; i < S.route.size(); i++) {
		tem_capacity = 0;
		for (auto j = 0; j < S.route[i].size(); j++) {
			tem_capacity += demand[S.route[i][j]];
		}
		tem_extra_capacity += max(0.0, tem_capacity - load_capacity);
	}
	if (abs(tem_extra_capacity - F.extra_capacity) > 0.01) {
		isExtraCapacity = false;
		std::cout << "△ Warning! The EXTRA CAPACITY is wrongly calculated! Original: " << F.extra_capacity << " Validation: " << tem_extra_capacity << endl;
	}

	
	bool isCumulativeTime = true;
	double temp_cumulative_time = 0;
	for (auto i = 0; i < S.route.size(); i++) {
		for (auto j = 1; j < S.route[i].size() - 1; j++) {
			temp_cumulative_time += S.arrival_time[i][j];
		}
	}
	if (abs(temp_cumulative_time - F.cumulative_time) > 0.01) {
		isCumulativeTime = false;
		std::cout << "△ Warning! The CUMULATIVE TIME is wrongly calculated! Original: " << F.cumulative_time << " Validation: " << temp_cumulative_time << endl;
	}


	
	if (isCapacity && isExtraCapacity && isCumulativeTime) {
		std::cout << endl;
		std::cout << "Congratulations!! All validations are passed!!" << endl;
	}

}


void deltaEvaluation(tempInformation& single, evaluationFunction& new_F, const int i) {
	arrivalTimeForSingleRoute(single.single_arrival_time, single.single_route);
	cumulativeTimeForSingleRoute(single.single_arrival_time, single.single_route, single.cumulative_time);
	double delta_cumulative_time = single.cumulative_time - S.single_cumulative_time[i];

	new_F.cumulative_time = F.cumulative_time + delta_cumulative_time;
	//new_F.extra_capacity = F.extra_capacity;
	new_F.evaluation_obj = new_F.cumulative_time + penalty.alpha_extra_capacity * new_F.extra_capacity + penalty.alpha_extra_vehicle * new_F.extra_vehicle;
}


void deltaEvaluation(tempInformation& single1, evaluationFunction& new_F, const int i, tempInformation& single2) {
	arrivalTimeForSingleRoute(single1.single_arrival_time, single1.single_route);
	cumulativeTimeForSingleRoute(single1.single_arrival_time, single1.single_route, single1.cumulative_time);
	double delta_cumulative_time1 = single1.cumulative_time - S.single_cumulative_time[i];

	arrivalTimeForSingleRoute(single2.single_arrival_time, single2.single_route);
	cumulativeTimeForSingleRoute(single2.single_arrival_time, single2.single_route, single2.cumulative_time);
	double delta_cumulative_time2 = single2.cumulative_time;

	new_F.cumulative_time = F.cumulative_time + delta_cumulative_time1 + delta_cumulative_time2;
	new_F.evaluation_obj = new_F.cumulative_time + penalty.alpha_extra_capacity * new_F.extra_capacity + penalty.alpha_extra_vehicle * new_F.extra_vehicle;
}

void deltaEvaluationForTwoSingle(tempInformation& single1, tempInformation& single2, evaluationFunction& new_F, const int i, const int j) {
	arrivalTimeForSingleRoute(single1.single_arrival_time, single1.single_route);
	cumulativeTimeForSingleRoute(single1.single_arrival_time, single1.single_route, single1.cumulative_time);
	double delta_cumulative_time1 = single1.cumulative_time - S.single_cumulative_time[i];


	arrivalTimeForSingleRoute(single2.single_arrival_time, single2.single_route);
	cumulativeTimeForSingleRoute(single2.single_arrival_time, single2.single_route, single2.cumulative_time);
	double delta_cumulative_time2 = single2.cumulative_time - S.single_cumulative_time[j];

	new_F.cumulative_time = F.cumulative_time + delta_cumulative_time1 + delta_cumulative_time2;
	new_F.evaluation_obj = new_F.cumulative_time + penalty.alpha_extra_capacity * new_F.extra_capacity + penalty.alpha_extra_vehicle * new_F.extra_vehicle;
}


// /* ================ operator ================ */


void relocate(routeInformation& S) {
	operatorPara para;
	para.select_route_index = randomInt(S.route.size());

	const auto& route = S.route[para.select_route_index];
	if (route.size() < 4) {
		//cout << "relocate random number is infeasible!" << endl;
		return;
	}

	para.select_node_index = randomInt(1, route.size() - 2);
	int node = route[para.select_node_index];

	if (Tabu[node][para.select_route_index] >= t) return;

	int new_index = randomInt(1, route.size() - 2);                  
	while (new_index == para.select_node_index) {
		new_index = randomInt(1, route.size() - 2);
	}

	
	double current_cum_time = S.single_cumulative_time[para.select_route_index];
	double new_cum_time = 0;
	double current_time = 0;
	int prev_node = 0;
	int old_pos = para.select_node_index;
	int new_pos = new_index;
	int u = node;

	
	int size = route.size();
	for (int k = 1; k < size; ++k) {
		int node_k;
		if (old_pos < new_pos) {
			
			if (k < old_pos) node_k = route[k];
			else if (k < new_pos) node_k = route[k + 1];
			else if (k == new_pos) node_k = u;
			else node_k = route[k];
		}
		else {
			
			if (k < new_pos) node_k = route[k];
			else if (k == new_pos) node_k = u;
			else if (k <= old_pos) node_k = route[k - 1];
			else node_k = route[k];
		}

		current_time += travelTime[prev_node][node_k];
		if (k < size - 1) {
			new_cum_time += current_time;
		}
		prev_node = node_k;
	}

	double delta_obj = new_cum_time - current_cum_time;
	

	if (F.evaluation_obj + delta_obj < F.evaluation_obj - 1e-9) {
		
		vector<int>& r = S.route[para.select_route_index];
		if (old_pos < new_pos) {
			rotate(r.begin() + old_pos, r.begin() + old_pos + 1, r.begin() + new_index + 1);
		}
		else {
			rotate(r.begin() + new_index, r.begin() + old_pos, r.begin() + old_pos + 1);
		}

		
		S.single_cumulative_time[para.select_route_index] = new_cum_time;
		update_route_arrival_time(S, para.select_route_index); 

		F.cumulative_time += delta_obj;
		F.evaluation_obj += delta_obj;

		
		Tabu[node][para.select_route_index] = tabu_step + t;
	}
}


void opt(routeInformation& S) {
	operatorPara para;
	para.select_route_index = randomInt(S.route.size());

	const auto& route = S.route[para.select_route_index];
	if (route.size() < 4) {
		//cout << "opt random number is infeasible!" << endl;
		return;
	}
	para.node1_index = randomInt(1, route.size() - 2);
	para.node2_index = randomInt(1, route.size() - 2);
	while (para.node1_index == para.node2_index) {
		para.node2_index = randomInt(1, route.size() - 2);
	}

	
	if (para.node1_index > para.node2_index) {
		int tem = para.node1_index;
		para.node1_index = para.node2_index;
		para.node2_index = tem;
	}

	
	double current_cum_time = S.single_cumulative_time[para.select_route_index];
	double new_cum_time = 0;
	double current_time = 0;
	int prev_node = 0;

	
	int size = route.size();
	for (int k = 1; k < size; ++k) {
		int node_k;
		if (k >= para.node1_index && k <= para.node2_index) {
			
			node_k = route[para.node2_index - (k - para.node1_index)];
		}
		else {
			node_k = route[k];
		}

		current_time += travelTime[prev_node][node_k];
		if (k < size - 1) {
			new_cum_time += current_time;
		}
		prev_node = node_k;
	}

	double delta_obj = new_cum_time - current_cum_time;

	if (F.evaluation_obj + delta_obj < F.evaluation_obj - 1e-9) {
		vector<int>& r = S.route[para.select_route_index];
		reverse(r.begin() + para.node1_index, r.begin() + para.node2_index + 1);

		S.single_cumulative_time[para.select_route_index] = new_cum_time;
		update_route_arrival_time(S, para.select_route_index);

		F.cumulative_time += delta_obj;
		F.evaluation_obj += delta_obj;
	}
}


void exchange(routeInformation& S) {
	int node1_index;
	int node2_index;
	operatorPara para;
	para.select_route_index = randomInt(S.route.size());

	const auto& route = S.route[para.select_route_index];
	if (route.size() < 4) {
		//cout << "exchange random number is infeasible!" << endl;
		return;
	}

	node1_index = randomInt(1, route.size() - 2);
	node2_index = randomInt(1, route.size() - 2);
	while (node1_index == node2_index) {
		node2_index = randomInt(1, route.size() - 2);
	}

	int node1 = route[node1_index];
	int node2 = route[node2_index];

	if (Tabu[node1][para.select_route_index] >= t || Tabu[node2][para.select_route_index] >= t) return;

	
	double current_cum_time = S.single_cumulative_time[para.select_route_index];
	double new_cum_time = 0;
	double current_time = 0;
	int prev_node = 0;

	
	int size = route.size();
	for (int k = 1; k < size; ++k) {
		int node_k;
		if (k == node1_index) node_k = node2;
		else if (k == node2_index) node_k = node1;
		else node_k = route[k];

		current_time += travelTime[prev_node][node_k];
		if (k < size - 1) {
			new_cum_time += current_time;
		}
		prev_node = node_k;
	}

	double delta_obj = new_cum_time - current_cum_time;

	if (F.evaluation_obj + delta_obj < F.evaluation_obj - 1e-9) {
		vector<int>& r = S.route[para.select_route_index];
		swap(r[node1_index], r[node2_index]);

		S.single_cumulative_time[para.select_route_index] = new_cum_time;
		update_route_arrival_time(S, para.select_route_index);

		F.cumulative_time += delta_obj;
		F.evaluation_obj += delta_obj;

		Tabu[node1][para.select_route_index] = tabu_step + t;
		Tabu[node2][para.select_route_index] = tabu_step + t;
	}
}


void relocate_sec(routeInformation& S) {
	operatorPara para;
	para.isChange = false;
	vector<int> available_route;
	for (auto i = 0; i < S.route.size(); i++) {
		if (S.route[i].size() - 2 >= 3) {
			available_route.push_back(i);
		}
	}

	if (available_route.size() > 0) {
		para.isChange = true;
		para.select_route_index = available_route[rand() % available_route.size()];
	}

	if (!para.isChange) return;

	const auto& route = S.route[para.select_route_index];
	para.select_node_index = randomInt(1, route.size() - 3);

	int last_node = route[para.select_node_index];
	int next_node = route[para.select_node_index + 1];

	if (Tabu[last_node][para.select_route_index] >= t || Tabu[next_node][para.select_route_index] >= t) return;

	// Random new position
	

	int len_after_remove = route.size() - 2;
	int new_index = randomInt(1, len_after_remove - 1);

	
	if (new_index == para.select_node_index) return;

	
	double current_cum_time = S.single_cumulative_time[para.select_route_index];
	double new_cum_time = 0;
	double current_time = 0;
	int prev_node = 0;
	int old_pos = para.select_node_index;
	int new_pos = new_index; // index in the route *after removal*

	

	// We iterate k from 0 to len_after_remove-1.
	// And insert u1, u2 at new_pos.

	for (int k = 1; k < len_after_remove; ++k) {
		// Determine node at step k of shortened route
		int node_k;
		if (k < old_pos) node_k = route[k];
		else node_k = route[k + 2];

		// Before processing node_k, check if we insert u1, u2 here
		if (k == new_pos) {
			current_time += travelTime[prev_node][last_node];
			new_cum_time += current_time;
			prev_node = last_node;

			current_time += travelTime[prev_node][next_node];
			new_cum_time += current_time;
			prev_node = next_node;
		}

		current_time += travelTime[prev_node][node_k];
		if (k < len_after_remove - 1) { // Not the last depot of shortened route
			new_cum_time += current_time;
		}
		prev_node = node_k;
	}

	double delta_obj = new_cum_time - current_cum_time;

	if (F.evaluation_obj + delta_obj < F.evaluation_obj - 1e-9) {
		vector<int>& r = S.route[para.select_route_index];
		vector<int> temp = { last_node, next_node };
		r.erase(r.begin() + old_pos, r.begin() + old_pos + 2);
		r.insert(r.begin() + new_pos, temp.begin(), temp.end());

		S.single_cumulative_time[para.select_route_index] = new_cum_time;
		update_route_arrival_time(S, para.select_route_index);

		F.cumulative_time += delta_obj;
		F.evaluation_obj += delta_obj;

		Tabu[last_node][para.select_route_index] = tabu_step + t;
		Tabu[next_node][para.select_route_index] = tabu_step + t;
	}
}


void swap_one_one(routeInformation& S) {
	operatorPara para;
	para.route1_index = randomInt(S.route.size());
	para.route2_index = randomInt(S.route.size());
	while (para.route1_index == para.route2_index) {
		para.route2_index = randomInt(S.route.size());
	}

	const auto& r1 = S.route[para.route1_index];
	const auto& r2 = S.route[para.route2_index];

	if (r1.size() < 3 || r2.size() < 3) {
		//cout << "swap_one_one random number is infeasible!" << endl;
		return;
	}

	para.node1_index = randomInt(1, r1.size() - 2);
	para.node2_index = randomInt(1, r2.size() - 2);

	int node1 = r1[para.node1_index];
	int node2 = r2[para.node2_index];

	if (Tabu[node1][para.route2_index] >= t || Tabu[node2][para.route1_index] >= t) return;

	// O(1) Delta Evaluation
	// Route 1 change
	int p1 = r1[para.node1_index - 1];
	int n1 = r1[para.node1_index + 1];
	double t_p1 = S.arrival_time[para.route1_index][para.node1_index - 1];
	double t_node1_old = S.arrival_time[para.route1_index][para.node1_index];

	double t_node2_new_at_r1 = t_p1 + travelTime[p1][node2];
	double delta_time_at_pos1 = t_node2_new_at_r1 - t_node1_old;

	double delta_shift_1 = (t_node2_new_at_r1 + travelTime[node2][n1]) - (t_node1_old + travelTime[node1][n1]);

	int customers_after_1 = r1.size() - 2 - para.node1_index;
	double delta_cum_1 = delta_time_at_pos1 + delta_shift_1 * customers_after_1;

	// Route 2 change
	int p2 = r2[para.node2_index - 1];
	int n2 = r2[para.node2_index + 1];
	double t_p2 = S.arrival_time[para.route2_index][para.node2_index - 1];
	double t_node2_old = S.arrival_time[para.route2_index][para.node2_index];

	double t_node1_new_at_r2 = t_p2 + travelTime[p2][node1];
	double delta_time_at_pos2 = t_node1_new_at_r2 - t_node2_old;

	double delta_shift_2 = (t_node1_new_at_r2 + travelTime[node1][n2]) - (t_node2_old + travelTime[node2][n2]);

	int customers_after_2 = r2.size() - 2 - para.node2_index;
	double delta_cum_2 = delta_time_at_pos2 + delta_shift_2 * customers_after_2;

	// Capacity & Penalty
	double delta_cum_total = delta_cum_1 + delta_cum_2;

	int cap1_old = S.used_capacity[para.route1_index];
	int cap2_old = S.used_capacity[para.route2_index];
	int cap1_new = cap1_old - demand[node1] + demand[node2];
	int cap2_new = cap2_old - demand[node2] + demand[node1];

	double delta_penalty = calc_penalty_delta(cap1_old, cap1_new - cap1_old, load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(cap2_old, cap2_new - cap2_old, load_capacity, penalty.alpha_extra_capacity);

	double total_delta = delta_cum_total + delta_penalty;

	if (total_delta < -1e-9) {
		// Apply
		S.route[para.route1_index][para.node1_index] = node2;
		S.route[para.route2_index][para.node2_index] = node1;

		S.used_capacity[para.route1_index] = cap1_new;
		S.used_capacity[para.route2_index] = cap2_new;

		update_route_arrival_time(S, para.route1_index);
		update_route_arrival_time(S, para.route2_index);
		cumulativeTimeForSingleRoute(S, para.route1_index);
		cumulativeTimeForSingleRoute(S, para.route2_index);

		F.cumulative_time += delta_cum_total;
		F.extra_capacity += (max(0.0, (double)cap1_new - load_capacity) - max(0.0, (double)cap1_old - load_capacity)) +
			(max(0.0, (double)cap2_new - load_capacity) - max(0.0, (double)cap2_old - load_capacity));
		F.evaluation_obj += total_delta;

		Tabu[node1][para.route2_index] = tabu_step + t;
		Tabu[node2][para.route1_index] = tabu_step + t;
	}
}


void swap_one_two(routeInformation& S) {
	operatorPara para;
	para.isChange = false;
	vector<int> available_route;
	for (auto i = 0; i < S.route.size(); i++) {
		if (S.route[i].size() > 3) {
			available_route.push_back(i);
		}
	}

	if (available_route.size() > 1) {
		para.isChange = true;
		para.route1_index = available_route[rand() % available_route.size()];
		para.route2_index = available_route[rand() % available_route.size()];
		while (para.route2_index == para.route1_index) {
			para.route2_index = available_route[rand() % available_route.size()];
		}
	}

	if (!para.isChange) return;

	const auto& r1 = S.route[para.route1_index];
	const auto& r2 = S.route[para.route2_index];

	if (r1.size() - 2 < 1 || r2.size() - 3 < 1) {
		//cout << "swap_one_two random number is infeasible!" << endl;
		return;
	}

	para.node1_index = randomInt(1, r1.size() - 2);        
	para.node2_index = randomInt(1, r2.size() - 3);        
	auto node1 = r1[para.node1_index];
	auto last_node2 = r2[para.node2_index];
	auto next_node2 = r2[para.node2_index + 1];

	if (Tabu[node1][para.node2_index] >= t || Tabu[last_node2][para.route1_index] >= t || Tabu[next_node2][para.route1_index] >= t) return;

	// O(1) Delta Evaluation

	// Route 1: Remove node1, Insert last_node2, next_node2 at same pos
	int p1 = r1[para.node1_index - 1];
	int n1 = r1[para.node1_index + 1];
	double t_p1 = S.arrival_time[para.route1_index][para.node1_index - 1];
	double t_node1_old = S.arrival_time[para.route1_index][para.node1_index];

	// New nodes arrival
	double t_new_1a = t_p1 + travelTime[p1][last_node2];
	double t_new_1b = t_new_1a + travelTime[last_node2][next_node2];

	double shift1 = (t_new_1b + travelTime[next_node2][n1]) - (t_node1_old + travelTime[node1][n1]);
	int customers_after_1 = r1.size() - 2 - para.node1_index;

	double delta_cum_1 = (t_new_1a + t_new_1b) - t_node1_old + shift1 * customers_after_1;

	// Route 2: Remove last_node2, next_node2, Insert node1 at same pos
	int p2 = r2[para.node2_index - 1];
	int n2 = r2[para.node2_index + 2];
	double t_p2 = S.arrival_time[para.route2_index][para.node2_index - 1];
	double t_node2a_old = S.arrival_time[para.route2_index][para.node2_index];
	double t_node2b_old = S.arrival_time[para.route2_index][para.node2_index + 1];

	double t_new_2 = t_p2 + travelTime[p2][node1];

	double shift2 = (t_new_2 + travelTime[node1][n2]) - (t_node2b_old + travelTime[next_node2][n2]);
	int customers_after_2 = r2.size() - 2 - (para.node2_index + 1);

	double delta_cum_2 = t_new_2 - (t_node2a_old + t_node2b_old) + shift2 * customers_after_2;

	double delta_cum_total = delta_cum_1 + delta_cum_2;

	// Capacity
	int cap1_old = S.used_capacity[para.route1_index];
	int cap2_old = S.used_capacity[para.route2_index];
	int cap1_new = cap1_old - demand[node1] + demand[last_node2] + demand[next_node2];
	int cap2_new = cap2_old - demand[last_node2] - demand[next_node2] + demand[node1];

	double delta_penalty = calc_penalty_delta(cap1_old, cap1_new - cap1_old, load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(cap2_old, cap2_new - cap2_old, load_capacity, penalty.alpha_extra_capacity);

	double total_delta = delta_cum_total + delta_penalty;

	if (total_delta < -1e-9) {
		vector<int> temp_node2 = { last_node2, next_node2 };

		// Apply R1
		S.route[para.route1_index].erase(S.route[para.route1_index].begin() + para.node1_index);
		S.route[para.route1_index].insert(S.route[para.route1_index].begin() + para.node1_index, temp_node2.begin(), temp_node2.end());
		S.used_capacity[para.route1_index] = cap1_new;

		// Apply R2
		S.route[para.route2_index].erase(S.route[para.route2_index].begin() + para.node2_index, S.route[para.route2_index].begin() + para.node2_index + 2);
		S.route[para.route2_index].insert(S.route[para.route2_index].begin() + para.node2_index, node1);
		S.used_capacity[para.route2_index] = cap2_new;

		update_route_arrival_time(S, para.route1_index);
		update_route_arrival_time(S, para.route2_index);
		cumulativeTimeForSingleRoute(S, para.route1_index);
		cumulativeTimeForSingleRoute(S, para.route2_index);

		F.cumulative_time += delta_cum_total;
		F.extra_capacity += (max(0.0, (double)cap1_new - load_capacity) - max(0.0, (double)cap1_old - load_capacity)) +
			(max(0.0, (double)cap2_new - load_capacity) - max(0.0, (double)cap2_old - load_capacity));
		F.evaluation_obj += total_delta;

		Tabu[node1][para.route2_index] = tabu_step + t;
		Tabu[last_node2][para.route1_index] = tabu_step + t;
		Tabu[next_node2][para.route1_index] = tabu_step + t;
	}
}


void swap_two_two(routeInformation& S) {
	operatorPara para;
	para.isChange = false;
	vector<int> available_route;
	for (auto i = 0; i < S.route.size(); i++) {
		if (S.route[i].size() > 3) {
			available_route.push_back(i);
		}
	}

	if (available_route.size() > 1) {
		para.isChange = true;
		para.route1_index = available_route[rand() % available_route.size()];
		para.route2_index = available_route[rand() % available_route.size()];
		while (para.route2_index == para.route1_index) {
			para.route2_index = available_route[rand() % available_route.size()];
		}
	}

	if (!para.isChange) return;

	const auto& r1 = S.route[para.route1_index];
	const auto& r2 = S.route[para.route2_index];

	if (r1.size() - 3 < 1 || r2.size() - 3 < 1) {
		//cout << "swap_two_two random number is infeasible!" << endl;
		return;
	}

	para.node1_index = randomInt(1, r1.size() - 3);        
	para.node2_index = randomInt(1, r2.size() - 3);        
	auto last_node1 = r1[para.node1_index];
	auto next_node1 = r1[para.node1_index + 1];
	auto last_node2 = r2[para.node2_index];
	auto next_node2 = r2[para.node2_index + 1];

	if (Tabu[last_node1][para.node2_index] >= t || Tabu[next_node1][para.node2_index] >= t || Tabu[last_node2][para.route1_index] >= t || Tabu[next_node2][para.route1_index] >= t) return;

	// O(1) Delta Evaluation

	// Route 1: Swap (last_node1, next_node1) with (last_node2, next_node2)
	int p1 = r1[para.node1_index - 1];
	int n1 = r1[para.node1_index + 2];
	double t_p1 = S.arrival_time[para.route1_index][para.node1_index - 1];
	double t_u1_old = S.arrival_time[para.route1_index][para.node1_index];
	double t_u2_old = S.arrival_time[para.route1_index][para.node1_index + 1];

	double t_v1_new = t_p1 + travelTime[p1][last_node2];
	double t_v2_new = t_v1_new + travelTime[last_node2][next_node2];

	double shift1 = (t_v2_new + travelTime[next_node2][n1]) - (t_u2_old + travelTime[next_node1][n1]);
	int cust_after_1 = r1.size() - 1 - (para.node1_index + 2);

	double delta_cum_1 = (t_v1_new + t_v2_new) - (t_u1_old + t_u2_old) + shift1 * cust_after_1;

	// Route 2: Swap (last_node2, next_node2) with (last_node1, next_node1)
	int p2 = r2[para.node2_index - 1];
	int n2 = r2[para.node2_index + 2];
	double t_p2 = S.arrival_time[para.route2_index][para.node2_index - 1];
	double t_v1_old = S.arrival_time[para.route2_index][para.node2_index];
	double t_v2_old = S.arrival_time[para.route2_index][para.node2_index + 1];

	double t_u1_new = t_p2 + travelTime[p2][last_node1];
	double t_u2_new = t_u1_new + travelTime[last_node1][next_node1];

	double shift2 = (t_u2_new + travelTime[next_node1][n2]) - (t_v2_old + travelTime[next_node2][n2]);
	int cust_after_2 = r2.size() - 1 - (para.node2_index + 2);

	double delta_cum_2 = (t_u1_new + t_u2_new) - (t_v1_old + t_v2_old) + shift2 * cust_after_2;

	double delta_cum_total = delta_cum_1 + delta_cum_2;

	// Capacity
	int cap1_old = S.used_capacity[para.route1_index];
	int cap2_old = S.used_capacity[para.route2_index];
	int cap1_new = cap1_old - demand[last_node1] - demand[next_node1] + demand[last_node2] + demand[next_node2];
	int cap2_new = cap2_old - demand[last_node2] - demand[next_node2] + demand[last_node1] + demand[next_node1];

	double delta_penalty = calc_penalty_delta(cap1_old, cap1_new - cap1_old, load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(cap2_old, cap2_new - cap2_old, load_capacity, penalty.alpha_extra_capacity);

	double total_delta = delta_cum_total + delta_penalty;

	if (total_delta < -1e-9) {
		vector<int> temp_node1 = { last_node1, next_node1 };
		vector<int> temp_node2 = { last_node2, next_node2 };

		S.route[para.route1_index].erase(S.route[para.route1_index].begin() + para.node1_index, S.route[para.route1_index].begin() + para.node1_index + 2);
		S.route[para.route1_index].insert(S.route[para.route1_index].begin() + para.node1_index, temp_node2.begin(), temp_node2.end());
		S.used_capacity[para.route1_index] = cap1_new;

		S.route[para.route2_index].erase(S.route[para.route2_index].begin() + para.node2_index, S.route[para.route2_index].begin() + para.node2_index + 2);
		S.route[para.route2_index].insert(S.route[para.route2_index].begin() + para.node2_index, temp_node1.begin(), temp_node1.end());
		S.used_capacity[para.route2_index] = cap2_new;

		update_route_arrival_time(S, para.route1_index);
		update_route_arrival_time(S, para.route2_index);
		cumulativeTimeForSingleRoute(S, para.route1_index);
		cumulativeTimeForSingleRoute(S, para.route2_index);

		F.cumulative_time += delta_cum_total;
		F.extra_capacity += (max(0.0, (double)cap1_new - load_capacity) - max(0.0, (double)cap1_old - load_capacity)) +
			(max(0.0, (double)cap2_new - load_capacity) - max(0.0, (double)cap2_old - load_capacity));
		F.evaluation_obj += total_delta;

		Tabu[last_node1][para.route2_index] = tabu_step + t;
		Tabu[next_node1][para.route2_index] = tabu_step + t;
		Tabu[last_node2][para.route1_index] = tabu_step + t;
		Tabu[next_node2][para.route1_index] = tabu_step + t;
	}
}


void shift(routeInformation& S) {
	operatorPara para;
	para.select_route_index = randomInt(S.route.size());                  

	if (S.route[para.select_route_index].size() - 2 < 1) {
		//cout << "shift random number is infeasible!" << endl;
		return;
	}

	para.select_node_index = randomInt(1, S.route[para.select_route_index].size() - 2);    

	int insert_route_index = randomInt(S.route.size());              
	
	while ((S.route[para.select_route_index].size() == 3 && insert_route_index == S.route.size()) || (insert_route_index == para.select_route_index)) {
		insert_route_index = randomInt(S.route.size());
	}

	int node = S.route[para.select_route_index][para.select_node_index];

	if (Tabu[node][insert_route_index] >= t) { return; }                 

	
	pair<double, double> delta_rem = calc_remove_delta(S, para.select_route_index, para.select_node_index);

	pair<double, double> delta_ins;
	int insert_pos = 0;
	bool is_new_route = (insert_route_index == S.route.size());

	if (is_new_route) {
		
		double arrival = travelTime[0][node];
		delta_ins = { (double)demand[node], arrival };
	}
	else {
		insert_pos = randomInt(1, S.route[insert_route_index].size() - 1);
		delta_ins = calc_insert_delta(S, insert_route_index, insert_pos, node);
	}

	double delta_cum_total = delta_rem.second + delta_ins.second;

	// Capacity
	int cap1_old = S.used_capacity[para.select_route_index];
	int cap1_new = cap1_old - demand[node];

	int cap2_old = is_new_route ? 0 : S.used_capacity[insert_route_index];
	int cap2_new = cap2_old + demand[node];

	double delta_penalty = calc_penalty_delta(cap1_old, cap1_new - cap1_old, load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(cap2_old, cap2_new - cap2_old, load_capacity, penalty.alpha_extra_capacity);

	// Vehicle Count
	double delta_veh_penalty = 0;
	if (is_new_route) {
		if (S.route.size() + 1 > vehicle_num) delta_veh_penalty += penalty.alpha_extra_vehicle;
	}
	else {
		
		if (S.route[para.select_route_index].size() == 3) {
			if (S.route.size() > vehicle_num) delta_veh_penalty -= penalty.alpha_extra_vehicle;
		}
	}

	double total_delta = delta_cum_total + delta_penalty + delta_veh_penalty;

	if (total_delta < -1e-9) {
		// Apply

		// 1. Insert into r2
		if (is_new_route) {
			S.route.push_back({ 0, node, 0 });
			S.used_capacity.push_back(demand[node]);
			S.arrival_time.emplace_back();
			S.single_cumulative_time.push_back(0); // updated below
			insert_route_index = S.route.size() - 1;
		}
		else {
			S.route[insert_route_index].insert(S.route[insert_route_index].begin() + insert_pos, node);
			S.used_capacity[insert_route_index] = cap2_new;
		}

		// 2. Remove from r1
		
		S.route[para.select_route_index].erase(S.route[para.select_route_index].begin() + para.select_node_index);
		S.used_capacity[para.select_route_index] = cap1_new;

		// 3. Clean up empty route
		bool r1_removed = false;
		if (S.route[para.select_route_index].size() <= 2) {
			S.route.erase(S.route.begin() + para.select_route_index);
			S.used_capacity.erase(S.used_capacity.begin() + para.select_route_index);
			S.arrival_time.erase(S.arrival_time.begin() + para.select_route_index);
			S.single_cumulative_time.erase(S.single_cumulative_time.begin() + para.select_route_index);
			r1_removed = true;

			
			if (insert_route_index > para.select_route_index) {
				insert_route_index--;
			}
		}
		else {
			update_route_arrival_time(S, para.select_route_index);
			cumulativeTimeForSingleRoute(S, para.select_route_index);
		}

		update_route_arrival_time(S, insert_route_index);
		cumulativeTimeForSingleRoute(S, insert_route_index);

		// Update Global F
		F.cumulative_time += delta_cum_total;
		// Re-calculate Extra Capacity/Vehicle to be safe (or maintain delta)
		// delta maintenance is faster but error-prone with vector erase.
		// Given erase is rare (only when size==3), re-calc F parts is ok.
		// But let's try to maintain delta.
		F.extra_capacity += (max(0.0, (double)cap1_new - load_capacity) - max(0.0, (double)cap1_old - load_capacity)) +
			(max(0.0, (double)cap2_new - load_capacity) - max(0.0, (double)cap2_old - load_capacity));

		int current_veh = S.route.size();
		F.extra_vehicle = max(0, current_veh - vehicle_num);

		F.evaluation_obj += total_delta; // This might drift slightly, but verify() will catch it.
		// Better: F.evaluation_obj = F.cumulative_time + penalty...

		Tabu[node][para.select_route_index] = tabu_step + t;
		Tabu[node][insert_route_index] = tabu_step + t;
	}
}


void shift_sec(routeInformation& S) {
	operatorPara para;
	para.isChange = false;
	vector<int> available_route;
	for (size_t i = 0; i < S.route.size(); i++) {
		if (S.route[i].size() > 3) {
			available_route.push_back(static_cast<int>(i));
		}
	}

	if (available_route.size() >= 1) {
		para.isChange = true;
		para.select_route_index = available_route[rand() % available_route.size()];
	}

	if (!para.isChange) { return; }

	if (S.route[para.select_route_index].size() - 3 < 1) {
		//cout << "shift_sec random number is infeasible!" << endl;
		return;
	}

	para.select_node_index = randomInt(1, S.route[para.select_route_index].size() - 3);

	//int insert_route_index = randomInt(S.route.size() + 1);              
	int insert_route_index = randomInt(S.route.size());              
	
	while ((S.route[para.select_route_index].size() == 4 && insert_route_index == S.route.size()) || (insert_route_index == para.select_route_index)) {
		//insert_route_index = randomInt(S.route.size() + 1);
		insert_route_index = randomInt(S.route.size());
	}

	auto last_node = S.route[para.select_route_index][para.select_node_index];
	auto next_node = S.route[para.select_route_index][para.select_node_index + 1];

	if (Tabu[last_node][insert_route_index] >= t || Tabu[next_node][insert_route_index] >= t) { return; }                

	

	// 1. Remove from R1
	int r1_idx = para.select_route_index;
	int pos1 = para.select_node_index;
	const auto& r1 = S.route[r1_idx];
	int p1 = r1[pos1 - 1];
	int n1 = r1[pos1 + 2];
	double t_last = S.arrival_time[r1_idx][pos1];
	double t_next = S.arrival_time[r1_idx][pos1 + 1];

	double shift_rem = travelTime[p1][n1] - (travelTime[p1][last_node] + travelTime[last_node][next_node] + travelTime[next_node][n1]);
	int nodes_after_rem = r1.size() - 1 - (pos1 + 2); // Customers after the pair

	double delta_rem = -t_last - t_next + shift_rem * nodes_after_rem;

	// 2. Insert into R2
	double delta_ins = 0;
	bool is_new_route = (insert_route_index == S.route.size());
	int insert_pos = 0;

	if (is_new_route) {
		double t1 = travelTime[0][last_node];
		double t2 = t1 + travelTime[last_node][next_node];
		delta_ins = t1 + t2; // + travelTime[next_node][0] doesn't count for cumulative
	}
	else {
		insert_pos = randomInt(1, S.route[insert_route_index].size() - 1);
		const auto& r2 = S.route[insert_route_index];
		int p2 = r2[insert_pos - 1];
		int n2 = r2[insert_pos];

		double t_p2 = S.arrival_time[insert_route_index][insert_pos - 1];
		double t1 = t_p2 + travelTime[p2][last_node];
		double t2 = t1 + travelTime[last_node][next_node];

		double shift_ins = (travelTime[p2][last_node] + travelTime[last_node][next_node] + travelTime[next_node][n2]) - travelTime[p2][n2];
		int nodes_after_ins = r2.size() - 1 - insert_pos;

		delta_ins = t1 + t2 + shift_ins * nodes_after_ins;
	}

	double delta_cum_total = delta_rem + delta_ins;

	// Capacity
	int cap1_old = S.used_capacity[r1_idx];
	int cap1_new = cap1_old - demand[last_node] - demand[next_node];

	int cap2_old = is_new_route ? 0 : S.used_capacity[insert_route_index];
	int cap2_new = cap2_old + demand[last_node] + demand[next_node];

	double delta_penalty = calc_penalty_delta(cap1_old, cap1_new - cap1_old, load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(cap2_old, cap2_new - cap2_old, load_capacity, penalty.alpha_extra_capacity);

	// Vehicle Count
	double delta_veh_penalty = 0;
	if (is_new_route) {
		if (S.route.size() + 1 > vehicle_num) delta_veh_penalty += penalty.alpha_extra_vehicle;
	}
	else {
		
		if (r1.size() == 4) {
			if (S.route.size() > vehicle_num) delta_veh_penalty -= penalty.alpha_extra_vehicle;
		}
	}

	double total_delta = delta_cum_total + delta_penalty + delta_veh_penalty;

	if (total_delta < -1e-9) {
		// Apply
		vector<int> temp_node = { last_node, next_node };

		// 1. Insert
		if (is_new_route) {
			S.route.push_back({ 0, last_node, next_node, 0 });
			S.used_capacity.push_back(demand[last_node] + demand[next_node]);
			S.arrival_time.emplace_back();
			S.single_cumulative_time.push_back(0);
			insert_route_index = S.route.size() - 1;
		}
		else {
			S.route[insert_route_index].insert(S.route[insert_route_index].begin() + insert_pos, temp_node.begin(), temp_node.end());
			S.used_capacity[insert_route_index] = cap2_new;
		}

		// 2. Remove
		S.route[r1_idx].erase(S.route[r1_idx].begin() + pos1, S.route[r1_idx].begin() + pos1 + 2);
		S.used_capacity[r1_idx] = cap1_new;

		// 3. Cleanup
		if (S.route[r1_idx].size() <= 2) {
			S.route.erase(S.route.begin() + r1_idx);
			S.used_capacity.erase(S.used_capacity.begin() + r1_idx);
			S.arrival_time.erase(S.arrival_time.begin() + r1_idx);
			S.single_cumulative_time.erase(S.single_cumulative_time.begin() + r1_idx);

			if (insert_route_index > r1_idx) insert_route_index--;
		}
		else {
			update_route_arrival_time(S, r1_idx);
			cumulativeTimeForSingleRoute(S, r1_idx);
		}

		update_route_arrival_time(S, insert_route_index);
		cumulativeTimeForSingleRoute(S, insert_route_index);

		F.cumulative_time += delta_cum_total;
		F.extra_capacity += (max(0, cap1_new - load_capacity) - max(0, cap1_old - load_capacity)) +
			(max(0, cap2_new - load_capacity) - max(0, cap2_old - load_capacity));
		int current_veh = S.route.size();
		F.extra_vehicle = max(0, current_veh - vehicle_num);
		F.evaluation_obj += total_delta;

		Tabu[last_node][para.select_route_index] = tabu_step + t;
		Tabu[last_node][insert_route_index] = tabu_step + t;
		Tabu[next_node][para.select_route_index] = tabu_step + t;
		Tabu[next_node][insert_route_index] = tabu_step + t;
	}
}


void opt_star(routeInformation& S) {
	operatorPara para;
	para.route1_index = randomInt(S.route.size());
	para.route2_index = randomInt(S.route.size());
	while (para.route1_index == para.route2_index) {
		para.route2_index = randomInt(S.route.size());
	}

	const auto& r1 = S.route[para.route1_index];
	const auto& r2 = S.route[para.route2_index];

	if (r1.size() < 3 || r2.size() < 3) {
		//cout << "opt_star random number is infeasible!" << endl;
		return;
	}

	para.node1_index = randomInt(1, r1.size() - 2);
	para.node2_index = randomInt(1, r2.size() - 2);

	if (para.node1_index == 1 && para.node2_index == 1) { return; }

	// O(N) Evaluation without vector copy

	// New R1: r1[0..node1_index] + r2[node2_index..end]
	// Actually opt* usually swaps tails: r1[0..i] connects to r2[j+1..end].
	// Current code: 
	// segment1: r1.begin + node1_index to end. (So removing node1_index itself? No, begin+idx includes idx)
	// Original code: vector<int> segment1(single1.single_route.begin() + para.node1_index, single1.single_route.end());
	// So tail starts AT node1_index.
	// New R1: r1[0 .. node1_index-1] + r2[node2_index .. end].

	// Let's verify original logic:
	// erase(begin + node1_index, end) -> keeps 0 .. node1_index-1.
	// insert(end, segment2) -> appends r2[node2_index .. end].

	// So New R1 = r1[0..p1-1] + r2[p2..end]
	// New R2 = r2[0..p2-1] + r1[p1..end]

	double new_cum_1 = 0;
	double cur_time_1 = 0;
	int prev_1 = 0;

	int cap1_new = 0;

	// Calc New R1
	// Part 1: r1 prefix
	for (int k = 0; k < para.node1_index; ++k) {
		int node = r1[k];
		if (k > 0) {
			cur_time_1 += travelTime[prev_1][node];
			new_cum_1 += cur_time_1;
			cap1_new += demand[node];
		}
		prev_1 = node;
	}
	// Part 2: r2 suffix
	for (int k = para.node2_index; k < r2.size(); ++k) {
		int node = r2[k];
		cur_time_1 += travelTime[prev_1][node];
		if (k < r2.size() - 1) { // Not last depot
			new_cum_1 += cur_time_1;
			cap1_new += demand[node];
		}
		prev_1 = node;
	}

	double new_cum_2 = 0;
	double cur_time_2 = 0;
	int prev_2 = 0;
	int cap2_new = 0;

	// Calc New R2
	// Part 1: r2 prefix
	for (int k = 0; k < para.node2_index; ++k) {
		int node = r2[k];
		if (k > 0) {
			cur_time_2 += travelTime[prev_2][node];
			new_cum_2 += cur_time_2;
			cap2_new += demand[node];
		}
		prev_2 = node;
	}
	// Part 2: r1 suffix
	for (int k = para.node1_index; k < r1.size(); ++k) {
		int node = r1[k];
		cur_time_2 += travelTime[prev_2][node];
		if (k < r1.size() - 1) {
			new_cum_2 += cur_time_2;
			cap2_new += demand[node];
		}
		prev_2 = node;
	}

	double delta_cum = (new_cum_1 + new_cum_2) - (S.single_cumulative_time[para.route1_index] + S.single_cumulative_time[para.route2_index]);

	int cap1_old = S.used_capacity[para.route1_index];
	int cap2_old = S.used_capacity[para.route2_index];

	double delta_penalty = calc_penalty_delta(cap1_old, cap1_new - cap1_old, load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(cap2_old, cap2_new - cap2_old, load_capacity, penalty.alpha_extra_capacity);

	double total_delta = delta_cum + delta_penalty;

	if (total_delta < -1e-9) {
		// Apply
		// Use temp vectors to swap safely
		vector<int> tail1(r1.begin() + para.node1_index, r1.end());
		vector<int> tail2(r2.begin() + para.node2_index, r2.end());

		S.route[para.route1_index].erase(S.route[para.route1_index].begin() + para.node1_index, S.route[para.route1_index].end());
		S.route[para.route1_index].insert(S.route[para.route1_index].end(), tail2.begin(), tail2.end());

		S.route[para.route2_index].erase(S.route[para.route2_index].begin() + para.node2_index, S.route[para.route2_index].end());
		S.route[para.route2_index].insert(S.route[para.route2_index].end(), tail1.begin(), tail1.end());

		S.used_capacity[para.route1_index] = cap1_new;
		S.used_capacity[para.route2_index] = cap2_new;

		S.single_cumulative_time[para.route1_index] = new_cum_1;
		S.single_cumulative_time[para.route2_index] = new_cum_2;

		update_route_arrival_time(S, para.route1_index);
		update_route_arrival_time(S, para.route2_index);

		F.cumulative_time += delta_cum;
		F.extra_capacity += (max(0, cap1_new - load_capacity) - max(0, cap1_old - load_capacity)) +
			(max(0, cap2_new - load_capacity) - max(0, cap2_old - load_capacity));
		F.evaluation_obj += total_delta;
	}
}


void dEAX(const routeInformation& S1, const routeInformation& S2, routeInformation& child, evaluationFunction& childF) {
	
	vector<vector<int>> nextA(node_num + 1);
	vector<vector<int>> nextB(node_num + 1);

	for (const auto& r : S1.route) {
		if (r.size() < 2) continue;
		for (size_t i = 0; i < r.size() - 1; ++i) {
			nextA[r[i]].push_back(r[i + 1]);
		}
	}
	for (const auto& r : S2.route) {
		if (r.size() < 2) continue;
		for (size_t i = 0; i < r.size() - 1; ++i) {
			nextB[r[i]].push_back(r[i + 1]);
		}
	}

	
	while (nextA[0].size() < nextB[0].size()) {
		nextA[0].push_back(0);
	}
	while (nextB[0].size() < nextA[0].size()) {
		nextB[0].push_back(0);
	}

	
	for (int u = 0; u <= node_num; ++u) {
		for (size_t i = 0; i < nextA[u].size(); ) {
			int v = nextA[u][i];
			auto itB = find(nextB[u].begin(), nextB[u].end(), v);
			if (itB != nextB[u].end()) {
				nextA[u].erase(nextA[u].begin() + i);
				nextB[u].erase(itB);
			}
			else {
				++i;
			}
		}
	}

	
	vector<vector<int>> prevB(node_num + 1);
	for (int u = 0; u <= node_num; ++u) {
		for (int v : nextB[u]) {
			prevB[v].push_back(u);
		}
	}

	
	struct Cycle {
		vector<pair<int, int>> A_edges;
		vector<pair<int, int>> B_edges;
		vector<int> nodes;
	};
	vector<Cycle> ab_cycles;

	for (int start = 0; start <= node_num; ++start) {
		while (!nextA[start].empty()) {
			Cycle cycle;
			int curr = start;
			while (true) {
				if (nextA[curr].empty()) break;
				int nxt = nextA[curr].back();
				nextA[curr].pop_back();

				cycle.A_edges.push_back({ curr, nxt });
				cycle.nodes.push_back(curr);
				cycle.nodes.push_back(nxt);

				if (prevB[nxt].empty()) break;
				int prev = prevB[nxt].back();
				prevB[nxt].pop_back();

				auto itB = find(nextB[prev].begin(), nextB[prev].end(), nxt);
				if (itB != nextB[prev].end()) nextB[prev].erase(itB);

				cycle.B_edges.push_back({ prev, nxt });
				cycle.nodes.push_back(prev);

				curr = prev;
				if (curr == start) break;
			}
			
			sort(cycle.nodes.begin(), cycle.nodes.end());
			cycle.nodes.erase(unique(cycle.nodes.begin(), cycle.nodes.end()), cycle.nodes.end());
			ab_cycles.push_back(cycle);
		}
	}

	
	if (ab_cycles.empty()) {
		child = S1;
		calInformation(child, childF);
		return;
	}

	int center_idx = randomInt(ab_cycles.size());
	vector<int> e_set_indices;
	e_set_indices.push_back(center_idx);

	vector<bool> center_has_node(node_num + 1, false);
	for (int n : ab_cycles[center_idx].nodes) center_has_node[n] = true;

	for (size_t i = 0; i < ab_cycles.size(); ++i) {
		if (i == center_idx) continue;
		bool share = false;
		for (int n : ab_cycles[i].nodes) {
			if (center_has_node[n]) {
				share = true;
				break;
			}
		}
		if (share) e_set_indices.push_back(i);
	}

	
	vector<vector<int>> adj(node_num + 1);
	for (const auto& r : S1.route) {
		if (r.size() < 2) continue;
		for (size_t i = 0; i < r.size() - 1; ++i) {
			adj[r[i]].push_back(r[i + 1]);
		}
	}

	
	while (adj[0].size() < nextA[0].size()) {
		adj[0].push_back(0);
	}

	for (int idx : e_set_indices) {
		for (const auto& edge : ab_cycles[idx].A_edges) {
			auto it = find(adj[edge.first].begin(), adj[edge.first].end(), edge.second);
			if (it != adj[edge.first].end()) adj[edge.first].erase(it);
		}
		for (const auto& edge : ab_cycles[idx].B_edges) {
			adj[edge.first].push_back(edge.second);
		}
	}

	
	vector<vector<int>> main_routes;
	vector<vector<int>> subtours;
	vector<bool> visited(node_num + 1, false);

	
	while (!adj[0].empty()) {
		vector<int> route;
		route.push_back(0);
		int curr = 0;
		while (true) {
			if (adj[curr].empty()) break;
			int nxt = adj[curr].back();
			adj[curr].pop_back();
			route.push_back(nxt);
			visited[nxt] = true;
			curr = nxt;
			if (curr == 0) break;
		}
		main_routes.push_back(route);
	}

	
	for (int i = 1; i <= node_num; ++i) {
		if (!visited[i] && !adj[i].empty()) {
			vector<int> subtour;
			int curr = i;
			while (true) {
				subtour.push_back(curr);
				visited[curr] = true;
				if (adj[curr].empty()) break;
				int nxt = adj[curr].back();
				adj[curr].pop_back();
				curr = nxt;
				if (curr == i) break;
			}
			subtours.push_back(subtour);
		}
	}

	
	for (const auto& subtour : subtours) {
		if (subtour.empty()) continue;

		int k = subtour.size();
		int subtour_cap = 0;
		for (int node : subtour) subtour_cap += demand[node];

		vector<double> SumInternal(k, 0.0);
		vector<double> PathTime(k, 0.0);
		double T_cycle_total = 0;
		for (int i = 0; i < k; ++i) {
			T_cycle_total += travelTime[subtour[i]][subtour[(i + 1) % k]];
		}

		for (int i = 0; i < k; ++i) { // i is the index of u, v is (i+1)%k
			int u = subtour[i];
			int v = subtour[(i + 1) % k];
			PathTime[i] = T_cycle_total - travelTime[u][v];

			double current_path_time = 0;
			double sum_internal = 0;
			int curr_idx = (i + 1) % k;
			for (int step = 1; step < k; ++step) {
				int nxt_idx = (curr_idx + 1) % k;
				current_path_time += travelTime[subtour[curr_idx]][subtour[nxt_idx]];
				sum_internal += current_path_time;
				curr_idx = nxt_idx;
			}
			SumInternal[i] = sum_internal;
		}

		double best_cost_diff = 1e18;
		int best_route_idx = -1;
		int best_sub_edge_idx = -1;
		int best_route_edge_idx = -1;

		for (size_t r = 0; r < main_routes.size(); ++r) {
			const auto& target_route = main_routes[r];
			int m = target_route.size();
			if (m < 2) continue;

			
			vector<double> arr_time(m, 0.0);
			double current_time = 0;
			int current_cap = 0;
			for (int j = 0; j < m - 1; ++j) {
				if (target_route[j] != 0) current_cap += demand[target_route[j]];
				current_time += travelTime[target_route[j]][target_route[j + 1]];
				arr_time[j + 1] = current_time;
			}

			double current_penalty = max(0.0, (double)(current_cap - load_capacity)) * penalty.alpha_extra_capacity;
			double new_penalty = max(0.0, (double)(current_cap + subtour_cap - load_capacity)) * penalty.alpha_extra_capacity;
			double delta_penalty = new_penalty - current_penalty;

			for (int j = 0; j < m - 1; ++j) {
				int x = target_route[j];
				int y = target_route[j + 1];
				double T_prev = arr_time[j];
				int N_after = m - 2 - j;
				if (N_after < 0) N_after = 0;

				for (int i = 0; i < k; ++i) {
					int u = subtour[i];
					int v = subtour[(i + 1) % k];

					double Delay = travelTime[x][v] + PathTime[i] + travelTime[u][y] - travelTime[x][y];
					
					double delta_CCVRP = k * (T_prev + travelTime[x][v]) + SumInternal[i] + Delay * N_after;

					double total_delta = delta_CCVRP + delta_penalty;

					if (total_delta < best_cost_diff) {
						best_cost_diff = total_delta;
						best_route_idx = r;
						best_sub_edge_idx = i;
						best_route_edge_idx = j;
					}
				}
			}
		}

		
		if (best_route_idx != -1 && best_sub_edge_idx != -1 && best_route_edge_idx != -1) {
			auto& target_route = main_routes[best_route_idx];
			vector<int> to_insert;
			to_insert.reserve(k);
			for (int step = 0; step < k; ++step) {
				to_insert.push_back(subtour[(best_sub_edge_idx + 1 + step) % k]);
			}
			target_route.insert(target_route.begin() + best_route_edge_idx + 1, to_insert.begin(), to_insert.end());
		}
	}

	
	vector<vector<int>> final_routes;
	for (const auto& r : main_routes) {
		if (r.size() > 2) { 
			final_routes.push_back(r);
		}
	}

	
	while (final_routes.size() < vehicle_num) {
		final_routes.push_back({ 0, 0 });
	}

	child.route = final_routes;
	calInformation(child, childF);
}


void swap_star(routeInformation& S) {
	operatorPara para;
	para.route1_index = randomInt(S.route.size());

	// Neighborhood guided selection
	int r1_idx = para.route1_index;
	int r2_idx = -1;

	// 50% chance: Random (Diversification)
	// 50% chance: Neighbor-guided (Intensification)
	if (randomDouble() < 0.5) {
		para.route2_index = randomInt(S.route.size());
		while (para.route1_index == para.route2_index) {
			para.route2_index = randomInt(S.route.size());
		}
	}
	else {
		if (S.route[r1_idx].size() > 2) {
			int node_idx = randomInt(1, S.route[r1_idx].size() - 2);
			int node = S.route[r1_idx][node_idx];

			if (sorted_neighbors.empty()) init_neighbors();

			// Try top neighbors
			for (int k = 0; k < min((int)sorted_neighbors[node].size(), 10); ++k) {
				int neighbor = sorted_neighbors[node][k];
				for (int r = 0; r < S.route.size(); ++r) {
					if (r == r1_idx) continue;
					for (int x : S.route[r]) {
						if (x == neighbor) {
							r2_idx = r;
							break;
						}
					}
					if (r2_idx != -1) break;
				}
				if (r2_idx != -1) break;
			}
		}

		if (r2_idx != -1) {
			para.route2_index = r2_idx;
		}
		else {
			para.route2_index = randomInt(S.route.size());
			while (para.route1_index == para.route2_index) {
				para.route2_index = randomInt(S.route.size());
			}
		}
	}

	const auto& r1 = S.route[para.route1_index];
	const auto& r2 = S.route[para.route2_index];

	if (r1.size() < 3 || r2.size() < 3) return;

	para.node1_index = randomInt(1, r1.size() - 2);
	para.node2_index = randomInt(1, r2.size() - 2);

	int node1 = r1[para.node1_index];
	int node2 = r2[para.node2_index];

	int max_routes = node_num * 2;
	if (para.route1_index >= max_routes || para.route2_index >= max_routes) return;
	if (Tabu[node1][para.route2_index] >= t || Tabu[node2][para.route1_index] >= t) return;

	
	auto delta_rem1 = calc_remove_delta(S, para.route1_index, para.node1_index);
	auto delta_rem2 = calc_remove_delta(S, para.route2_index, para.node2_index);

	double t_rem_cum1 = S.single_cumulative_time[para.route1_index] + delta_rem1.second;
	double t_rem_cum2 = S.single_cumulative_time[para.route2_index] + delta_rem2.second;

	int cap_rem1 = S.used_capacity[para.route1_index] - demand[node1];
	int cap_rem2 = S.used_capacity[para.route2_index] - demand[node2];

	
	double best_ins_cum2 = 1e18;
	int best_ins_pos2 = -1;
	int route2_len_after_rem = r2.size() - 1; 

	
	vector<int> virtual_r2 = r2;
	virtual_r2.erase(virtual_r2.begin() + para.node2_index);
	vector<double> virtual_arr2(route2_len_after_rem, 0.0);
	double current_time = 0;
	for (int i = 0; i < route2_len_after_rem - 1; ++i) {
		current_time += travelTime[virtual_r2[i]][virtual_r2[i + 1]];
		virtual_arr2[i + 1] = current_time;
	}

	for (int p = 1; p < route2_len_after_rem; ++p) {
		int prev = virtual_r2[p - 1];
		int next = virtual_r2[p];
		double t_prev = virtual_arr2[p - 1];

		double delta_dist = travelTime[prev][node1] + travelTime[node1][next] - travelTime[prev][next];
		double t_node1 = t_prev + travelTime[prev][node1];

		int nodes_after = route2_len_after_rem - 1 - p; 
		double delta_cum = t_node1 + delta_dist * nodes_after;

		if (delta_cum < best_ins_cum2) {
			best_ins_cum2 = delta_cum;
			best_ins_pos2 = p;
		}
	}

	
	double best_ins_cum1 = 1e18;
	int best_ins_pos1 = -1;
	int route1_len_after_rem = r1.size() - 1;

	vector<int> virtual_r1 = r1;
	virtual_r1.erase(virtual_r1.begin() + para.node1_index);
	vector<double> virtual_arr1(route1_len_after_rem, 0.0);
	current_time = 0;
	for (int i = 0; i < route1_len_after_rem - 1; ++i) {
		current_time += travelTime[virtual_r1[i]][virtual_r1[i + 1]];
		virtual_arr1[i + 1] = current_time;
	}

	for (int p = 1; p < route1_len_after_rem; ++p) {
		int prev = virtual_r1[p - 1];
		int next = virtual_r1[p];
		double t_prev = virtual_arr1[p - 1];

		double delta_dist = travelTime[prev][node2] + travelTime[node2][next] - travelTime[prev][next];
		double t_node2 = t_prev + travelTime[prev][node2];

		int nodes_after = route1_len_after_rem - 1 - p;
		double delta_cum = t_node2 + delta_dist * nodes_after;

		if (delta_cum < best_ins_cum1) {
			best_ins_cum1 = delta_cum;
			best_ins_pos1 = p;
		}
	}

	
	int cap1_new = cap_rem1 + demand[node2];
	int cap2_new = cap_rem2 + demand[node1];

	double delta_penalty = calc_penalty_delta(S.used_capacity[para.route1_index], cap1_new - S.used_capacity[para.route1_index], load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(S.used_capacity[para.route2_index], cap2_new - S.used_capacity[para.route2_index], load_capacity, penalty.alpha_extra_capacity);

	double old_cum = S.single_cumulative_time[para.route1_index] + S.single_cumulative_time[para.route2_index];
	double new_cum = t_rem_cum1 + best_ins_cum1 + t_rem_cum2 + best_ins_cum2;

	if (new_cum - old_cum + delta_penalty < -1e-9) {
		
		virtual_r1.insert(virtual_r1.begin() + best_ins_pos1, node2);
		virtual_r2.insert(virtual_r2.begin() + best_ins_pos2, node1);

		S.route[para.route1_index] = virtual_r1;
		S.route[para.route2_index] = virtual_r2;

		int cap1_old = S.used_capacity[para.route1_index];
		int cap2_old = S.used_capacity[para.route2_index];

		S.used_capacity[para.route1_index] = cap1_new;
		S.used_capacity[para.route2_index] = cap2_new;

		
		double old_single_cum1 = S.single_cumulative_time[para.route1_index];
		double old_single_cum2 = S.single_cumulative_time[para.route2_index];

		update_route_arrival_time(S, para.route1_index);
		update_route_arrival_time(S, para.route2_index);
		cumulativeTimeForSingleRoute(S, para.route1_index);
		cumulativeTimeForSingleRoute(S, para.route2_index);

		
		double real_delta_cum = (S.single_cumulative_time[para.route1_index] - old_single_cum1) +
			(S.single_cumulative_time[para.route2_index] - old_single_cum2);

		F.cumulative_time += real_delta_cum;
		F.extra_capacity += (max(0.0, (double)cap1_new - load_capacity) - max(0.0, (double)cap1_old - load_capacity)) +
			(max(0.0, (double)cap2_new - load_capacity) - max(0.0, (double)cap2_old - load_capacity));

		
		F.evaluation_obj = F.cumulative_time +
			penalty.alpha_extra_capacity * F.extra_capacity +
			penalty.alpha_extra_vehicle * F.extra_vehicle;

		Tabu[node1][para.route2_index] = tabu_step + t;
		Tabu[node2][para.route1_index] = tabu_step + t;
	}
}


void node_arc_exchange(routeInformation& S) {
	operatorPara para;
	para.route1_index = randomInt(S.route.size());
	para.route2_index = randomInt(S.route.size()); 

	const auto& r1 = S.route[para.route1_index];
	const auto& r2 = S.route[para.route2_index];

	if (r1.size() < 3 || r2.size() < 4) return;

	para.node1_index = randomInt(1, r1.size() - 2);
	para.node2_index = randomInt(1, r2.size() - 3); 

	
	if (para.route1_index == para.route2_index) {
		if (abs(para.node1_index - para.node2_index) <= 1 || para.node1_index == para.node2_index + 1) return;
	}

	int node = r1[para.node1_index];
	int arc1 = r2[para.node2_index];
	int arc2 = r2[para.node2_index + 1];

	int max_routes = node_num * 2;
	if (para.route1_index >= max_routes || para.route2_index >= max_routes) return;
	if (Tabu[node][para.route2_index] >= t || Tabu[arc1][para.route1_index] >= t || Tabu[arc2][para.route1_index] >= t) return;

	tempInformation single1, single2;
	evaluationFunction new_F = F;

	if (para.route1_index != para.route2_index) {
		single1.single_route = r1;
		single2.single_route = r2;

		vector<int> arc = { arc1, arc2 };
		single1.single_route.erase(single1.single_route.begin() + para.node1_index);
		single2.single_route.erase(single2.single_route.begin() + para.node2_index, single2.single_route.begin() + para.node2_index + 2);

		single1.single_route.insert(single1.single_route.begin() + para.node1_index, arc.begin(), arc.end());
		single2.single_route.insert(single2.single_route.begin() + para.node2_index, node);

		int cap1_new = S.used_capacity[para.route1_index] - demand[node] + demand[arc1] + demand[arc2];
		int cap2_new = S.used_capacity[para.route2_index] - demand[arc1] - demand[arc2] + demand[node];

		int cap1_old = S.used_capacity[para.route1_index];
		int cap2_old = S.used_capacity[para.route2_index];

		double delta_extra_cap = (max(0.0, (double)cap1_new - load_capacity) - max(0.0, (double)cap1_old - load_capacity)) +
			(max(0.0, (double)cap2_new - load_capacity) - max(0.0, (double)cap2_old - load_capacity));
		new_F.extra_capacity += delta_extra_cap;

		deltaEvaluationForTwoSingle(single1, single2, new_F, para.route1_index, para.route2_index);

		
		if (new_F.evaluation_obj < F.evaluation_obj - 1e-9) {
			S.route[para.route1_index] = single1.single_route;
			S.arrival_time[para.route1_index] = single1.single_arrival_time;
			S.single_cumulative_time[para.route1_index] = single1.cumulative_time;
			S.used_capacity[para.route1_index] = cap1_new;

			S.route[para.route2_index] = single2.single_route;
			S.arrival_time[para.route2_index] = single2.single_arrival_time;
			S.single_cumulative_time[para.route2_index] = single2.cumulative_time;
			S.used_capacity[para.route2_index] = cap2_new;

			F = new_F;
			Tabu[node][para.route2_index] = tabu_step + t;
			Tabu[arc1][para.route1_index] = tabu_step + t;
			Tabu[arc2][para.route1_index] = tabu_step + t;
		}
	}
	else {
		single1.single_route = r1;
		vector<int>& route = single1.single_route;

		int p1 = para.node1_index;
		int p2 = para.node2_index;

		int node_val = route[p1];
		int arc1_val = route[p2];
		int arc2_val = route[p2 + 1];

		
		vector<int> new_route;
		new_route.reserve(route.size());

		if (p1 < p2) {
			for (int i = 0; i < p1; ++i) new_route.push_back(route[i]);
			new_route.push_back(arc1_val);
			new_route.push_back(arc2_val);
			for (int i = p1 + 1; i < p2; ++i) new_route.push_back(route[i]);
			new_route.push_back(node_val);
			for (int i = p2 + 2; i < route.size(); ++i) new_route.push_back(route[i]);
		}
		else {
			for (int i = 0; i < p2; ++i) new_route.push_back(route[i]);
			new_route.push_back(node_val);
			for (int i = p2 + 2; i < p1; ++i) new_route.push_back(route[i]);
			new_route.push_back(arc1_val);
			new_route.push_back(arc2_val);
			for (int i = p1 + 1; i < route.size(); ++i) new_route.push_back(route[i]);
		}

		single1.single_route = new_route;

		deltaEvaluation(single1, new_F, para.route1_index);

		if (new_F.evaluation_obj < F.evaluation_obj - 1e-9) {
			S.route[para.route1_index] = single1.single_route;
			S.arrival_time[para.route1_index] = single1.single_arrival_time;
			S.single_cumulative_time[para.route1_index] = single1.cumulative_time;
			F = new_F;
		}
	}
}

// Ejection Chain Mutation 
void ejection_chain_mutation(routeInformation& S) {
	
	vector<int> valid_routes;
	for (int i = 0; i < S.route.size(); ++i) {
		if (S.route[i].size() >= 3) {
			valid_routes.push_back(i);
		}
	}

	if (valid_routes.size() < 3) return;

	
	int r1_idx = randomInt(valid_routes.size());
	int r2_idx = randomInt(valid_routes.size());
	while (r2_idx == r1_idx) r2_idx = randomInt(valid_routes.size());
	int r3_idx = randomInt(valid_routes.size());
	while (r3_idx == r1_idx || r3_idx == r2_idx) r3_idx = randomInt(valid_routes.size());

	int r1 = valid_routes[r1_idx];
	int r2 = valid_routes[r2_idx];
	int r3 = valid_routes[r3_idx];

	
	int p1 = randomInt(1, S.route[r1].size() - 2);
	int p2 = randomInt(1, S.route[r2].size() - 2);
	int p3 = randomInt(1, S.route[r3].size() - 2);

	int a = S.route[r1][p1];
	int b = S.route[r2][p2];
	int c = S.route[r3][p3];

	
	int max_routes = node_num * 2;
	if (r1 >= max_routes || r2 >= max_routes || r3 >= max_routes) return;
	if (Tabu[a][r2] >= t || Tabu[b][r3] >= t || Tabu[c][r1] >= t) return;

	
	int cap1_new = S.used_capacity[r1] - demand[a] + demand[c];
	int cap2_new = S.used_capacity[r2] - demand[b] + demand[a];
	int cap3_new = S.used_capacity[r3] - demand[c] + demand[b];

	double delta_penalty = calc_penalty_delta(S.used_capacity[r1], cap1_new - S.used_capacity[r1], load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(S.used_capacity[r2], cap2_new - S.used_capacity[r2], load_capacity, penalty.alpha_extra_capacity) +
		calc_penalty_delta(S.used_capacity[r3], cap3_new - S.used_capacity[r3], load_capacity, penalty.alpha_extra_capacity);

	
	double t_prev1 = S.arrival_time[r1][p1 - 1];
	double t_a_old = S.arrival_time[r1][p1];
	int prev1 = S.route[r1][p1 - 1];
	int next1 = S.route[r1][p1 + 1];
	double t_c_new = t_prev1 + travelTime[prev1][c];
	double shift1 = (t_c_new + travelTime[c][next1]) - (t_a_old + travelTime[a][next1]);
	int after1 = S.route[r1].size() - 2 - p1;
	double delta_cum1 = t_c_new - t_a_old + shift1 * after1;

	
	double t_prev2 = S.arrival_time[r2][p2 - 1];
	double t_b_old = S.arrival_time[r2][p2];
	int prev2 = S.route[r2][p2 - 1];
	int next2 = S.route[r2][p2 + 1];
	double t_a_new = t_prev2 + travelTime[prev2][a];
	double shift2 = (t_a_new + travelTime[a][next2]) - (t_b_old + travelTime[b][next2]);
	int after2 = S.route[r2].size() - 2 - p2;
	double delta_cum2 = t_a_new - t_b_old + shift2 * after2;

	
	double t_prev3 = S.arrival_time[r3][p3 - 1];
	double t_c_old = S.arrival_time[r3][p3];
	int prev3 = S.route[r3][p3 - 1];
	int next3 = S.route[r3][p3 + 1];
	double t_b_new = t_prev3 + travelTime[prev3][b];
	double shift3 = (t_b_new + travelTime[b][next3]) - (t_c_old + travelTime[c][next3]);
	int after3 = S.route[r3].size() - 2 - p3;
	double delta_cum3 = t_b_new - t_c_old + shift3 * after3;

	double total_delta = delta_cum1 + delta_cum2 + delta_cum3 + delta_penalty;

	
	if (total_delta < -1e-9) {
		S.route[r1][p1] = c;
		S.route[r2][p2] = a;
		S.route[r3][p3] = b;

		int cap1_old = S.used_capacity[r1];
		int cap2_old = S.used_capacity[r2];
		int cap3_old = S.used_capacity[r3];

		S.used_capacity[r1] = cap1_new;
		S.used_capacity[r2] = cap2_new;
		S.used_capacity[r3] = cap3_new;

		update_route_arrival_time(S, r1);
		update_route_arrival_time(S, r2);
		update_route_arrival_time(S, r3);

		
		double old_single_cum1 = S.single_cumulative_time[r1];
		double old_single_cum2 = S.single_cumulative_time[r2];
		double old_single_cum3 = S.single_cumulative_time[r3];

		cumulativeTimeForSingleRoute(S, r1);
		cumulativeTimeForSingleRoute(S, r2);
		cumulativeTimeForSingleRoute(S, r3);

		
		double real_delta_cum = (S.single_cumulative_time[r1] - old_single_cum1) +
			(S.single_cumulative_time[r2] - old_single_cum2) +
			(S.single_cumulative_time[r3] - old_single_cum3);

		F.cumulative_time += real_delta_cum;
		F.extra_capacity += (max(0.0, (double)cap1_new - load_capacity) - max(0.0, (double)cap1_old - load_capacity)) +
			(max(0.0, (double)cap2_new - load_capacity) - max(0.0, (double)cap2_old - load_capacity)) +
			(max(0.0, (double)cap3_new - load_capacity) - max(0.0, (double)cap3_old - load_capacity));

		
		F.evaluation_obj = F.cumulative_time +
			penalty.alpha_extra_capacity * F.extra_capacity +
			penalty.alpha_extra_vehicle * F.extra_vehicle;

		Tabu[c][r1] = tabu_step + t;
		Tabu[a][r2] = tabu_step + t;
		Tabu[b][r3] = tabu_step + t;
	}
}

// Route Reversal Operator (AVNS)
void route_reversal(routeInformation& S) {
	bool improvement = true;
	while (improvement) {
		improvement = false;
		for (int r = 0; r < S.route.size(); ++r) {
			if (S.route[r].size() <= 3) continue;

			double current_cost = S.single_cumulative_time[r];

			double reversed_cost = 0;
			double current_time = 0;
			int prev_node = 0;

			// Reverse: 0 -> vn -> vn-1 -> ... -> v1 -> 0
			// Iterate backwards from end-1 (vn) to 1 (v1)
			for (int k = S.route[r].size() - 2; k >= 1; --k) {
				int node = S.route[r][k];
				current_time += travelTime[prev_node][node];
				reversed_cost += current_time;
				prev_node = node;
			}
			// Add return to depot time (not added to cumulative)
			// current_time += travelTime[prev_node][0];

			if (reversed_cost < current_cost - 1e-9) {
				// Reverse only the customers (from index 1 to size-2)
				reverse(S.route[r].begin() + 1, S.route[r].end() - 1);

				// Update all info
				S.single_cumulative_time[r] = reversed_cost;
				update_route_arrival_time(S, r);

				// IMPORTANT: We must update F correctly.
				// Since F.cumulative_time is the sum of all single_cumulative_time,
				// and reversed_cost is smaller than current_cost:
				F.cumulative_time -= (current_cost - reversed_cost);

				// F.evaluation_obj must be re-calculated completely to avoid precision drift
				F.evaluation_obj = F.cumulative_time +
					penalty.alpha_extra_capacity * F.extra_capacity +
					penalty.alpha_extra_vehicle * F.extra_vehicle;

				improvement = true;
			}
		}
	}
}
static vector<int> index_NL = { 11, 12, 13, 14, 8, 9, 1, 3, 4, 2, 5, 6, 7 };


void RVNS() {

	FUNC_TIMER();          

	int r_NL, index_r_NL;
	int count = 0;
	int l_index_NL = index_NL.size();

	
	static vector<double> success_rate(15, 0.5); 
	static vector<int> success_count(15, 0);
	static vector<int> total_count(15, 1);

	
	static int shuffle_count = 0;
	shuffle_count++;
	
	if (shuffle_count % 10 == 0) {
		
		int mid_point = l_index_NL / 2;
		static std::mt19937 g(std::random_device{}());
		shuffle(index_NL.begin() + mid_point, index_NL.end(), g);
	}

	while (count < l_index_NL) {
		r_NL = index_NL[count];
		evaluationFunction old_F = F;
		switch (r_NL) {
		case 1:
			relocate(S); // 1-insertion
			break;
		case 2:
			exchange(S); // 1-1 exchange
			break;
		case 3:
			swap_one_one(S); // 1-1 swap
			break;
		case 4:
			shift(S); // 2-insertion
			break;
		case 5:
			relocate_sec(S); 
			break;
		case 6:
			swap_one_two(S); // 1-2 swap
			break;
		case 7:
			swap_two_two(S); // 2-2 swap
			break;
		case 8: {
			// 2-opt
			int opt_iterations = success_rate[8] > 0.3 ? 2 : 1;
			for (int i = 0; i < opt_iterations; i++) {
				opt(S);
			}
			break;
		}
		case 9: {
			// 2-opt*
			int opt_star_iterations = success_rate[9] > 0.3 ? 2 : 1;
			for (int i = 0; i < opt_star_iterations; i++) {
				opt_star(S);
			}
			break;
		}
		case 11: {
			int swap_star_iterations = success_rate[11] > 0.3 ? 2 : 1;
			for (int i = 0; i < swap_star_iterations; i++) {
				swap_star(S);
			}
			break;
		}
		case 12: {
			
			int n5_iterations = success_rate[12] > 0.3 ? 2 : 1;
			for (int i = 0; i < n5_iterations; i++) {
				node_arc_exchange(S);
			}
			break;
		}
		case 13: {
			
			int mut_iterations = success_rate[13] > 0.3 ? 2 : 1;
			for (int i = 0; i < mut_iterations; i++) {
				ejection_chain_mutation(S);
			}
			break;
		}
		case 14: {
			// Route Reversal (AVNS)
			route_reversal(S);
			break;
		}
		default:
			break;
		}

		total_count[r_NL]++;
		op_used_num[r_NL]++;
		if (F.evaluation_obj < old_F.evaluation_obj) {
			success_count[r_NL]++;
			op_effec_num[r_NL]++;
			
			count = 0;
		}
		else {
			count++;
		}

		
		success_rate[r_NL] = (double)success_count[r_NL] / total_count[r_NL];
	}

	
	static int adjust_count = 0;
	adjust_count++;
	if (adjust_count % 50 == 0) {
		
		vector<pair<double, int>> rate_op;
		for (int op : index_NL) {
			rate_op.emplace_back(success_rate[op], op);
		}
		
		sort(rate_op.rbegin(), rate_op.rend());
		
		int keep_ratio = (int)(index_NL.size() * 0.7);
		for (int i = 0; i < keep_ratio; i++) {
			index_NL[i] = rate_op[i].second;
		}
		
		vector<int> temp_ops;
		for (int i = keep_ratio; i < index_NL.size(); i++) {
			temp_ops.push_back(rate_op[i].second);
		}
		static std::mt19937 g(std::random_device{}());
		shuffle(temp_ops.begin(), temp_ops.end(), g);
		for (int i = 0; i < temp_ops.size(); i++) {
			index_NL[keep_ratio + i] = temp_ops[i];
		}
	}
}


// ALNS Parameters and Structures
struct ALNS_Parameters {
	double sigma1 = 50; // New global best
	double sigma2 = 20; // Better than current
	double sigma3 = 5;  // Accepted
	double reaction_factor = 0.1;

	vector<double> removal_weights;
	vector<double> insertion_weights;
	vector<double> removal_scores;
	vector<double> insertion_scores;
	vector<int> removal_counts;
	vector<int> insertion_counts;

	// Removals: 0: Random, 1: Worst, 2: Shaw(Time), 3: Shaw(Dist)
	// Insertions: 0: Greedy, 1: Regret-2
	int n_removal = 4;
	int n_insertion = 2;

	void init() {
		removal_weights.assign(n_removal, 1.0);
		insertion_weights.assign(n_insertion, 1.0);
		removal_scores.assign(n_removal, 0.0);
		insertion_scores.assign(n_insertion, 0.0);
		removal_counts.assign(n_removal, 0);
		insertion_counts.assign(n_insertion, 0);
	}

	// Roulette wheel selection
	int select_removal() {
		double total = 0;
		for (double w : removal_weights) total += w;
		double r = randomDouble() * total;
		double sum = 0;
		for (int i = 0; i < n_removal; i++) {
			sum += removal_weights[i];
			if (r <= sum) return i;
		}
		return n_removal - 1;
	}

	int select_insertion() {
		double total = 0;
		for (double w : insertion_weights) total += w;
		double r = randomDouble() * total;
		double sum = 0;
		for (int i = 0; i < n_insertion; i++) {
			sum += insertion_weights[i];
			if (r <= sum) return i;
		}
		return n_insertion - 1;
	}

	void update_scores(int r_op, int i_op, double score) {
		removal_scores[r_op] += score;
		insertion_scores[i_op] += score;
		removal_counts[r_op]++;
		insertion_counts[i_op]++;
	}

	void update_weights() {
		for (int i = 0; i < n_removal; i++) {
			if (removal_counts[i] > 0) {
				removal_weights[i] = removal_weights[i] * (1 - reaction_factor) + reaction_factor * (removal_scores[i] / removal_counts[i]);
				removal_scores[i] = 0;
				removal_counts[i] = 0;
			}
		}
		for (int i = 0; i < n_insertion; i++) {
			if (insertion_counts[i] > 0) {
				insertion_weights[i] = insertion_weights[i] * (1 - reaction_factor) + reaction_factor * (insertion_scores[i] / insertion_counts[i]);
				insertion_scores[i] = 0;
				insertion_counts[i] = 0;
			}
		}
	}
} alns_param;

bool alns_initialized = false;

// ALNS Helper Functions

// Calculate removal cost of a node (approximation for efficiency)
double calculate_removal_cost(const routeInformation& s, int route_idx, int node_idx) {
	// Cost reduction = current_obj - new_obj
	// Removing node i reduces obj by:
	// 1. Arrival time of i (t_i)
	// 2. Delay reduction for all subsequent nodes j: (t_i + t_{i->next}) - t_{prev->next}

	if (route_idx >= s.route.size() || node_idx <= 0 || node_idx >= s.route[route_idx].size() - 1) return 0;

	int node = s.route[route_idx][node_idx];
	int prev = s.route[route_idx][node_idx - 1];
	int next = s.route[route_idx][node_idx + 1];

	double t_i = s.arrival_time[route_idx][node_idx];
	double time_saved = t_i; // Node i is removed, so its arrival time is gone from sum

	// Calculate time shift for subsequent nodes
	double old_segment = travelTime[prev][node] + travelTime[node][next];


	double delta = (travelTime[prev][node] + travelTime[node][next]) - travelTime[prev][next];
	// Note: serviceTime[node] contributes to delay of next nodes.

	int nodes_after = s.route[route_idx].size() - node_idx - 2;
	/*int nodes_after = 0;
	for (int k = node_idx + 1; k < s.route[route_idx].size() - 1; k++) {
		nodes_after++;
	}*/
	return t_i + delta * nodes_after;
}


// vector<vector<int>> sorted_neighbors;

void init_neighbors() {
	if (!sorted_neighbors.empty()) return;
	sorted_neighbors.resize(customer_num + 1);

	for (int i = 1; i <= customer_num; i++) {
		vector<pair<double, int>> dists;
		for (int j = 1; j <= customer_num; j++) {
			if (i == j) continue;
			// Shaw relatedness: distance + |demand_diff|
			// Normalize? Or just raw sum.
			// Distance is usually large (e.g. 10-100), demand is 1-10.
			// Let's just use distance for general "closeness"
			double d = travelTime[i][j] + abs(demand[i] - demand[j]);
			dists.push_back({ d, j });
		}
		sort(dists.begin(), dists.end());
		for (auto& p : dists) {
			sorted_neighbors[i].push_back(p.second);
		}
	}
}


void apply_removal(routeInformation& s, const vector<bool>& is_removed) {
	// Rebuild routes skipping removed nodes
	// Also re-calculate arrival times and capacity
	// This is O(N) and much faster than repeated erase

	int write_r = 0;
	for (int r = 0; r < s.route.size(); r++) {
		vector<int>& route = s.route[r];
		vector<int> new_route;
		new_route.reserve(route.size());

		int current_load = 0;
		int prev = 0;

		// Always keep depot 0 at start
		new_route.push_back(0);

		bool has_customers = false;
		for (int i = 1; i < route.size() - 1; i++) {
			int node = route[i];
			if (!is_removed[node]) {
				new_route.push_back(node);
				current_load += demand[node];
				has_customers = true;
				prev = node;
			}
		}
		new_route.push_back(0); // End depot

		if (has_customers) {
			// Move to write position (compact routes)
			if (write_r != r) {
				s.route[write_r] = move(new_route);
			}
			else {
				s.route[write_r] = move(new_route); // self assignment optimized or same
			}
			s.used_capacity[write_r] = current_load;

			// Re-calc arrival times
			s.arrival_time[write_r].resize(s.route[write_r].size());
			s.arrival_time[write_r][0] = 0;
			double time = 0;
			int p_node = 0;
			for (int k = 1; k < s.route[write_r].size(); k++) {
				int curr = s.route[write_r][k];
				time += travelTime[p_node][curr];
				s.arrival_time[write_r][k] = time;
				p_node = curr;
			}

			// Recalc single cumulative
			double cum = 0;
			for (int k = 1; k < s.route[write_r].size() - 1; ++k) {
				cum += s.arrival_time[write_r][k];
			}
			s.single_cumulative_time[write_r] = cum;

			write_r++;
		}
	}

	// Resize to remove empty slots
	s.route.resize(write_r);
	s.used_capacity.resize(write_r);
	s.arrival_time.resize(write_r);
	s.single_cumulative_time.resize(write_r);
}

void remove_random(routeInformation& s, int q, vector<int>& removed, vector<bool>& is_removed) {
	// Randomly select q nodes that are not yet removed
	// To do this efficiently without retry loops:
	// Collect all available nodes
	static vector<int> available;
	available.clear();
	for (int r = 0; r < s.route.size(); r++) {
		for (int i = 1; i < s.route[r].size() - 1; i++) {
			int node = s.route[r][i];
			if (!is_removed[node]) {
				available.push_back(node);
			}
		}
	}

	if (available.empty()) return;

	static std::mt19937 g(std::random_device{}());
	shuffle(available.begin(), available.end(), g);

	int count = 0;
	for (int node : available) {
		if (count >= q) break;
		removed.push_back(node);
		is_removed[node] = true;
		count++;
	}
}

void remove_worst(routeInformation& s, int q, vector<int>& removed, vector<bool>& is_removed) {
	// Calculate removal cost for all NON-REMOVED nodes
	static vector<pair<double, int>> costs; // cost -> node_idx (global)
	costs.clear();

	// Need to know route and index for calc
	// Since we don't modify structure yet, indices are stable
	for (int r = 0; r < s.route.size(); r++) {
		for (int i = 1; i < s.route[r].size() - 1; i++) {
			int node = s.route[r][i];
			if (is_removed[node]) continue;

			double cost = calculate_removal_cost(s, r, i);
			double noise = randomDouble() * 0.2 + 0.9;
			costs.push_back({ cost * noise, node });
		}
	}

	sort(costs.rbegin(), costs.rend()); // High cost first

	for (int k = 0; k < min((int)costs.size(), q); k++) {
		int node = costs[k].second;
		removed.push_back(node);
		is_removed[node] = true;
	}
}

void remove_shaw(routeInformation& s, int q, vector<int>& removed, vector<bool>& is_removed) {
	if (sorted_neighbors.empty()) init_neighbors();

	// Pick random seed from available
	static vector<int> available;
	available.clear();
	for (int r = 0; r < s.route.size(); r++) {
		for (int i = 1; i < s.route[r].size() - 1; i++) {
			int node = s.route[r][i];
			if (!is_removed[node]) available.push_back(node);
		}
	}

	if (available.empty()) return;

	int seed = available[randomInt(available.size())];
	removed.push_back(seed);
	is_removed[seed] = true;

	while (removed.size() < q) {
		int r_node = removed[randomInt(removed.size())];

		// Find nearest available neighbor
		int best_node = -1;
		for (int neighbor : sorted_neighbors[r_node]) {
			if (!is_removed[neighbor]) {
				best_node = neighbor;
				break; // Found best
			}
		}

		if (best_node != -1) {
			removed.push_back(best_node);
			is_removed[best_node] = true;
		}
		else {
			// No neighbors found (all removed?), pick random
			// Re-scan available or just pick one if any left
			// For simplicity, just break or pick from remaining
			// Let's try to pick another random seed from remaining
			bool found = false;
			for (int i = 1; i <= customer_num; i++) {
				if (!is_removed[i]) {
					removed.push_back(i);
					is_removed[i] = true;
					found = true;
					break;
				}
			}
			if (!found) break;
		}
	}
}

// Optimized calculate_insertion_cost with pre-calculated nodes_after
double calculate_insertion_cost(const routeInformation& s, int route_idx, int pos_idx, int node, int nodes_after) {
	// Calculate cost increase
	if (route_idx >= s.route.size()) return 1e9;

	int prev = s.route[route_idx][pos_idx - 1];
	int next = s.route[route_idx][pos_idx];

	double delta = travelTime[prev][node] + travelTime[node][next] - travelTime[prev][next];

	// Arrival time of inserted node
	double t_prev = s.arrival_time[route_idx][pos_idx - 1];
	double t_arrive_node = t_prev + travelTime[prev][node]; // Approx
	double cost_increase = t_arrive_node + delta * nodes_after;

	return cost_increase;
}

// Wrapper for backward compatibility (if needed, but we should update calls)
double calculate_insertion_cost(const routeInformation& s, int route_idx, int pos_idx, int node) {
	int nodes_after = 0;
	for (int k = pos_idx; k < s.route[route_idx].size(); k++) {
		if (s.route[route_idx][k] != 0) nodes_after++;
	}
	return calculate_insertion_cost(s, route_idx, pos_idx, node, nodes_after);
}

// Helper to update arrival times for a single route
void update_route_arrival_time(routeInformation& S, int r) {
	int size = S.route[r].size();
	S.arrival_time[r].resize(size);
	double current_time = 0;
	S.arrival_time[r][0] = 0;

	for (int j = 0; j < size - 1; j++) {
		int u = S.route[r][j];
		int v = S.route[r][j + 1];
		current_time += travelTime[u][v];
		S.arrival_time[r][j + 1] = current_time;
	}
}


struct BestInsertInfo {
	int route_idx;
	int pos_idx;
	double cost;
};


void update_insertion_cache(const routeInformation& S, const vector<int>& removed_nodes,
	vector<bool>& node_inserted,
	vector<vector<pair<double, int>>>& cache_cost_pos, // [node_idx_in_removed][route_idx] -> {cost, pos}
	int specific_route_idx = -1)
{
	int start_r = 0;
	int end_r = S.route.size();
	if (specific_route_idx != -1) {
		start_r = specific_route_idx;
		end_r = specific_route_idx + 1;
	}

	for (int r = start_r; r < end_r; ++r) {
		
		int total_customers = 0;
		for (int x : S.route[r]) {
			if (x != 0) total_customers++;
		}
		int capacity_used = S.used_capacity[r];

		for (int i = 0; i < removed_nodes.size(); ++i) {
			if (node_inserted[i]) continue;

			int node = removed_nodes[i];
			if (capacity_used + demand[node] > load_capacity) {
				cache_cost_pos[i][r] = { 1e18, -1 };
				continue;
			}

			double best_cost = 1e18;
			int best_pos = -1;

			int customers_passed = 0;
			for (int p = 1; p < S.route[r].size(); ++p) {
				int prev = S.route[r][p - 1];
				int next = S.route[r][p];
				int nodes_after = total_customers - customers_passed;
				double t_prev = S.arrival_time[r][p - 1];

				
				double delta_dist = travelTime[prev][node] + travelTime[node][next] - travelTime[prev][next];
				double arrival_at_node = t_prev + travelTime[prev][node];
				double cost = arrival_at_node + delta_dist * nodes_after;

				if (cost < best_cost) {
					best_cost = cost;
					best_pos = p;
				}

				if (next != 0) customers_passed++;
			}
			cache_cost_pos[i][r] = { best_cost, best_pos };
		}
	}
}

void insert_greedy(routeInformation& s, vector<int>& removed) {
	if (removed.empty()) return;

	static vector<bool> node_inserted;
	node_inserted.assign(removed.size(), false);
	static vector<vector<pair<double, int>>> cache;

	if (cache.size() < removed.size()) {
		cache.resize(removed.size());
	}
	for (int i = 0; i < removed.size(); ++i) {
		if (cache[i].size() < s.route.size()) {
			cache[i].resize(s.route.size());
		}
	}

	
	update_insertion_cache(s, removed, node_inserted, cache, -1);

	for (int iter = 0; iter < removed.size(); ++iter) {
		double best_global_cost = 1e18;
		int best_node_idx = -1;
		int best_route_idx = -1;
		int best_pos_idx = -1;

		
		for (int i = 0; i < removed.size(); ++i) {
			if (node_inserted[i]) continue;

			
			for (int r = 0; r < s.route.size(); ++r) {
				if (cache[i][r].first < best_global_cost) {
					best_global_cost = cache[i][r].first;
					best_node_idx = i;
					best_route_idx = r;
					best_pos_idx = cache[i][r].second;
				}
			}

			
			// Cost = arrival_time + 0 (no nodes after) = travelTime[0][node]
			double new_route_cost = travelTime[0][removed[i]];
			
			if (s.route.size() >= vehicle_num) new_route_cost += 1e9; 

			if (new_route_cost < best_global_cost) {
				best_global_cost = new_route_cost;
				best_node_idx = i;
				best_route_idx = -1; // New route
				best_pos_idx = 1;
			}
		}

		if (best_node_idx == -1) break; // Should not happen

		
		int node = removed[best_node_idx];
		node_inserted[best_node_idx] = true;

		if (best_route_idx == -1) {
			
			s.route.push_back({ 0, node, 0 });
			s.used_capacity.push_back(demand[node]);
			s.arrival_time.emplace_back();
			s.single_cumulative_time.push_back(0); 
			update_route_arrival_time(s, s.route.size() - 1);
			cumulativeTimeForSingleRoute(s, s.route.size() - 1);

			
			for (auto& row : cache) {
				row.emplace_back(1e18, -1);
			}
			
			update_insertion_cache(s, removed, node_inserted, cache, s.route.size() - 1);
		}
		else {
			
			s.route[best_route_idx].insert(s.route[best_route_idx].begin() + best_pos_idx, node);
			s.used_capacity[best_route_idx] += demand[node];
			update_route_arrival_time(s, best_route_idx);
			cumulativeTimeForSingleRoute(s, best_route_idx);

			
			update_insertion_cache(s, removed, node_inserted, cache, best_route_idx);
		}
	}

	
	removed.clear();
}

void insert_regret(routeInformation& s, vector<int>& removed, int k_regret) {
	if (removed.empty()) return;

	static vector<bool> node_inserted;
	node_inserted.assign(removed.size(), false);
	static vector<vector<pair<double, int>>> cache;
	if (cache.size() < removed.size()) {
		cache.resize(removed.size());
	}
	for (int i = 0; i < removed.size(); ++i) {
		if (cache[i].size() < s.route.size()) {
			cache[i].resize(s.route.size());
		}
	}

	update_insertion_cache(s, removed, node_inserted, cache, -1);

	int remaining_nodes = removed.size();
	while (remaining_nodes > 0) {
		double max_regret = -1e18;
		int best_node_idx = -1;
		int best_route_idx = -1;
		int best_pos_idx = -1;

		for (int i = 0; i < removed.size(); ++i) {
			if (node_inserted[i]) continue;

			double best_cost = 1e18;
			double second_best_cost = 1e18;
			int cur_best_r = -1;
			int cur_best_p = -1;

			for (int r = 0; r < s.route.size(); ++r) {
				if (cache[i][r].second != -1) {
					double c = cache[i][r].first;
					if (c < best_cost) {
						second_best_cost = best_cost;
						best_cost = c;
						cur_best_r = r;
						cur_best_p = cache[i][r].second;
					}
					else if (c < second_best_cost) {
						second_best_cost = c;
					}
				}
			}

			
			double new_route_cost = travelTime[0][removed[i]];
			if (s.route.size() >= vehicle_num) new_route_cost += 1e9;
			if (new_route_cost < best_cost) {
				second_best_cost = best_cost;
				best_cost = new_route_cost;
				cur_best_r = -1;
				cur_best_p = 1;
			}
			else if (new_route_cost < second_best_cost) {
				second_best_cost = new_route_cost;
			}

			double regret = 0;
			if (second_best_cost < 1e17) {
				regret = second_best_cost - best_cost;
			}

			if (regret > max_regret) {
				max_regret = regret;
				best_node_idx = i;
				best_route_idx = cur_best_r;
				best_pos_idx = cur_best_p;
			}
		}

		if (best_node_idx == -1) break; // Should not happen

		
		int node = removed[best_node_idx];
		node_inserted[best_node_idx] = true;
		remaining_nodes--;

		if (best_route_idx == -1) {
			
			s.route.push_back({ 0, node, 0 });
			s.used_capacity.push_back(demand[node]);
			s.arrival_time.emplace_back();
			s.single_cumulative_time.push_back(0);
			update_route_arrival_time(s, s.route.size() - 1);
			cumulativeTimeForSingleRoute(s, s.route.size() - 1);

			
			for (auto& row : cache) {
				row.emplace_back(1e18, -1);
			}
			update_insertion_cache(s, removed, node_inserted, cache, s.route.size() - 1);
		}
		else {
			s.route[best_route_idx].insert(s.route[best_route_idx].begin() + best_pos_idx, node);
			s.used_capacity[best_route_idx] += demand[node];
			update_route_arrival_time(s, best_route_idx);
			cumulativeTimeForSingleRoute(s, best_route_idx);

			// Lazy Update
			update_insertion_cache(s, removed, node_inserted, cache, best_route_idx);
		}
	}
	removed.clear();
}

// Adaptive Large Neighborhood Search
void ALNS() {

	FUNC_TIMER();

	if (!alns_initialized) {
		alns_param.init();
		alns_initialized = true;
	}

	// 1. Select Operators
	int r_op = alns_param.select_removal();
	int i_op = alns_param.select_insertion();

	// 2. Determine number of customers to remove (gamma)
	// Random between 10% and 30% of customers, or min 4
	int min_rem = max(4, (int)(0.1 * customer_num));
	int max_rem = max(5, (int)(0.3 * customer_num));
	int q = randomInt(min_rem, max_rem + 1);

	// 3. Apply Removal (Logical)
	static vector<int> removed_nodes;
	removed_nodes.clear();
	static vector<bool> is_removed;
	is_removed.assign(customer_num + 1, false);

	switch (r_op) {
	case 0: remove_random(S, q, removed_nodes, is_removed); break;
	case 1: remove_worst(S, q, removed_nodes, is_removed); break;
	case 2:
	case 3: remove_shaw(S, q, removed_nodes, is_removed); break;
	}

	// Backup before modification
	routeInformation s_backup = S;
	evaluationFunction f_backup = F;

	// Apply Physical Removal & Re-calc basics
	apply_removal(S, is_removed);

	// 4. Apply Insertion
	switch (i_op) {
	case 0: insert_greedy(S, removed_nodes); break;
	case 1: insert_regret(S, removed_nodes, 2); break;
	}

	// 5. Update F (Incremental from components)
	F.cumulative_time = 0;
	F.extra_capacity = 0;
	for (int r = 0; r < S.route.size(); ++r) {
		F.cumulative_time += S.single_cumulative_time[r];
		F.extra_capacity += max(0.0, (double)S.used_capacity[r] - load_capacity);
	}
	F.extra_vehicle = max(0, (int)S.route.size() - vehicle_num);
	F.evaluation_obj = F.cumulative_time +
		F.extra_capacity * penalty.alpha_extra_capacity +
		F.extra_vehicle * penalty.alpha_extra_vehicle;

	// 6. Evaluation and Acceptance
	double score = 0;
	bool accept = false;

	// Update Global Best (Only if feasible and better)
	if (F.extra_capacity == 0 && F.extra_vehicle == 0 && F.cumulative_time < optimal.obj) {
		score = alns_param.sigma1;
		optimal.route = S.route;
		optimal.obj = F.cumulative_time;
		accept = true;
	}
	else if (F.evaluation_obj < f_backup.evaluation_obj) { // Better than current
		score = alns_param.sigma2;
		accept = true;
	}
	else {
		// Simulated Annealing Acceptance
		double temp = 100.0 * pow(0.999, t); // Simple cooling
		double delta = F.evaluation_obj - f_backup.evaluation_obj;
		if (randomDouble() < exp(-delta / temp)) {
			score = alns_param.sigma3;
			accept = true;
		}
	}

	if (!accept) {
		S = s_backup;
		F = f_backup;
	}

	// 7. Update Scores
	alns_param.update_scores(r_op, i_op, score);

	// Periodically update weights (every 50 calls)
	static int alns_iter = 0;
	alns_iter++;
	if (alns_iter % 50 == 0) {
		alns_param.update_weights();
	}
}

// VNS function incorporating RVNS and ALNS
void VNS() {
	// RVNS for local intensification
	RVNS();
	// ALNS for large neighborhood diversification
	ALNS();
}



int main()
{
	int c = 1;
	int c_max = 1;
	vector<double> cumulative_t;
	vector<double> op_time;
	while (c <= c_max) {
		t = 0;
		//readInitial_CMT(inFile);
		readVRP(inFile);
		cout << "The file name is: " << inFile << endl;
		getVehicleNum(vehicle_num_file, inst);
		cout << "The vehicle num is: " << vehicle_num << endl;
		calDistance();
		routeInformation elit;
		evaluationFunction elitF;
		GRASP_initialization(elit);
		S = elit;
		int it = p.max_initial_LS;
		while (it) {
			RVNS();
			it--;
		}
		elit = S;
		elitF = F;

		routeInformation newS1;
		GRASP_initialization(newS1);
		S = newS1;
		it = p.max_initial_LS;
		while (it) {
			RVNS();
			it--;
		}

		verify(S);
		verify(S, F);
		vector<double> convergence;
		

		startTime = clock();
		endTime = startTime;
		evaluationFunction best_F = F;
		routeInformation best_S = S;
		int no_improve = 0;
		int shake_no_improve = 0;
		int cycle = 0;

		double bestTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
		double global_best_time = 0;
		while (bestTime < timeMax) {

			VNS();

			if (F.evaluation_obj < best_F.evaluation_obj && S.route.size() <= vehicle_num) {
				if (best_F.evaluation_obj - F.evaluation_obj > 0.009) {
					endTime = clock();
					bestTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
					global_best_time = bestTime;
				}
				best_S = S;
				best_F = F;
				std::cout << std::fixed << std::setprecision(2);
				std::cout << "Iteration: " << t << ", Current best cumulative waiting time: " << best_F.cumulative_time << ", Taking time: " << global_best_time << "s, Current best evaluation objective: " << best_F.evaluation_obj << endl;
				no_improve = 0;
				shake_no_improve = 0;

				
				p.omega = p.min_omega;                             
				p.shake_s = max(2.0, p.shake_s - 0.5);             
				p.max_length_RCL_percent = max(0.10, p.max_length_RCL_percent - 0.02); 
			}
			else {
				no_improve++;
				shake_no_improve++;

				
				if (no_improve % 500 == 0) {
					p.omega = min(p.max_omega, p.omega + 1);             
					p.shake_s = min(8.0, p.shake_s + 0.5);               
					p.max_length_RCL_percent = min(0.30, p.max_length_RCL_percent + 0.02); 
				}
			}
			if (shake_no_improve >= p.max_shake_no_improve) {
				p.omega = p.max_omega;
				shake_no_improve = 0;
				
				p.shake_s = min(10.0, p.shake_s + 1.0);
			}
			else if (no_improve >= p.max_no_improve_new_S) {
				
				routeInformation child1, child2, child;
				evaluationFunction child1F, child2F, childF;
				dEAX(S, elit, child1, child1F);
				dEAX(elit, S, child2, child2F);

				if (child1F.evaluation_obj < child2F.evaluation_obj) {
					child = child1;
					childF = child1F;
				}
				else {
					child = child2;
					childF = child2F;
				}

				if (F.evaluation_obj < elitF.evaluation_obj && S.route.size() <= vehicle_num) {
					elit = S;
					elitF = F;
				}

				S = child;
				F = childF;

				cycle++;
				if (cycle % p.crossover_no_improve_cycle == 0) {
					routeInformation newS;
					GRASP_initialization(newS);
					S = newS;
					it = p.max_initial_LS;
					while (it) {
						RVNS();
						it--;
					}
				}

				no_improve = 0;
			}
			endTime = clock();
			bestTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
			//convergence.push_back(best_F.cumulative_time);
			no_improve++;
			shake_no_improve++;
			//convergence.push_back(best_F.cumulative_time);
			t++;

		}

		/*std::ofstream convergence_file("convergence_B-n68-k9/" + inst + "LLM-FAELS.txt", std::ios::app);
		if (convergence_file) {
			for (auto i = 0; i < convergence.size(); i++) {
				convergence_file << convergence[i] << " ";
			}
			convergence_file << endl;
			convergence_file << endl;
		}*/

		FunctionTimer::PrintReport();  

		verify(best_S);
		verify(best_S, best_F);
		std::cout << endl;
		std::cout << "best_F.cumulative_time: " << best_F.cumulative_time << endl;
		std::cout << "best_F.evaluation_obj: " << best_F.evaluation_obj << endl;
		std::cout << "The number of used vehicles is: " << best_S.route.size() << endl;
		std::cout << "The best time for the optimal solution is: " << global_best_time << endl;
		std::cout << "best_F.extra_vehicle: " << best_F.extra_vehicle << endl;
		std::cout << "The total iteration is: " << t << endl;
		std::cout << endl;



		cumulative_t.push_back(best_F.cumulative_time);
		op_time.push_back(global_best_time);
		c++;

		std::cout << endl;
		std::cout << "The final route is: " << endl;
		for (auto i = 0; i < best_S.route.size(); i++) {
			for (auto j = 0; j < best_S.route[i].size(); j++) {
				std::cout << best_S.route[i][j] << " ";
			}
			std::cout << endl;
		}
		std::cout << endl;

		


		
		/*std::ofstream outfile(instName + "result/" + inst + "_route.txt", std::ios::app);
		if (outfile) {
			outfile << "The final route of EX" << c - 1 << " is: \n" << endl;
			for (auto i = 0; i < best_S.route.size(); i++) {
				for (auto j = 0; j < best_S.route[i].size(); j++) {
					outfile << best_S.route[i][j] << " ";
				}
				outfile << endl;
			}
			outfile << endl;
		}*/
	}

	cout << "===================== Final results print =====================" << endl;
	cout << "Final optimal cumulative time of each run: " << endl;
	for (auto x : cumulative_t) {
		cout << x << " ";
	}
	cout << endl;
	cout << "Best operational time of each run: " << endl;
	for (auto x : op_time) {
		cout << x << " ";
	}
	cout << endl;

	
	/*cout << "The number of used times of operators:" << endl;
	for (auto i = 0; i < op_used_num.size(); i++) {
		cout << op_used_num[i] << " ";
	}
	cout << endl;
	cout << "The number of effective used times of operators:" << endl;
	for (auto i = 0; i < op_effec_num.size(); i++) {
		cout << op_effec_num[i] << " ";
	}
	cout << endl;*/

	
	/*std::ofstream outfile_time(instName + "result/" + inst + "_obj_runtime.txt", std::ios::app);
	if (outfile_time) {
		outfile_time << std::fixed << std::setprecision(2);
		outfile_time << "Final optimal cumulative time of each run: " << endl;
		for (auto x : cumulative_t) {
			outfile_time << x << endl;
		}
		outfile_time << endl;
		outfile_time << "Best operational time of each run: " << endl;
		for (auto x : op_time) {
			outfile_time << x << endl;
		}
		outfile_time << endl;
	}*/

	
	/*std::ofstream outfile_op("op_information.txt", std::ios::app);
	if (outfile_op) {
		outfile_op << "The number of used times of operators:" << endl;
		for (auto i = 0; i < op_used_num.size(); i++) {
			outfile_op << op_used_num[i] << endl;
		}
		outfile_op << endl;

		outfile_op << "The number of effective used times of operators:" << endl;
		for (auto i = 0; i < op_effec_num.size(); i++) {
			outfile_op << op_effec_num[i] << endl;
		}
		outfile_op << endl;
	}*/


	/*std::cout << endl;
	std::cout << "The final route is: " << endl;
	for (auto i = 0; i < best_S.route.size(); i++) {
		for (auto j = 0; j < best_S.route[i].size(); j++) {
			std::cout << best_S.route[i][j] << " ";
		}
		std::cout << endl;
	}
	std::cout << endl;*/

	/*cout << "The situation of the convergence is: " << endl;
	for (auto i = 0; i < convergence.size(); i++) {
		cout << convergence[i] << " ";
	}
	cout << endl;*/
}


