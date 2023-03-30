//============================================================================
// Name        : nauty-thesis.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

/*
 - Remove edges from graph
 - Find orbits
 - Get distinct list of orbits
 - Find largest and check if it is a used node
 */
#include <iostream>
#include <chrono>
#include <Common.h>
#include <vector>
#include <mutex>
#include <algorithm>
#include <initializer_list>


#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/adjacency_list.hpp>

struct GlobalCallback : public IloCplex::Callback::Function
{
	int& currentBound;
	GlobalCallback(int& currentBound);
	void invoke(const IloCplex::Callback::Context& context);
};
GlobalCallback::GlobalCallback(int& currentBoundIn) : currentBound(currentBoundIn)
{
}
void GlobalCallback::invoke(const IloCplex::Callback::Context& context)
{
	IloEnv env = context.getEnv();
	//float globalBound = context.getIncumbentObjective();
	//IloNumVar var(env);
	float globalBound = 9999;

	IloCplex::CplexStatus status = context.getRelaxationStatus(0);

	if (status == CPX_STAT_OPTIMAL){
		globalBound = context.getRelaxationObjective();
	}

	if (currentBound > (globalBound - 0.1)) {
		context.abort();
	}

	//if (currentBound > (bestBound + 0.1)) {
	//	context.abort();
	//}

}

int vertex_edge(int N, int i, int j){
	int upper_lim = N - 1;
	int lower_lim = N - i;
	int diff = j - i - 1;

	int edge = (i * (upper_lim + lower_lim)) / 2 + diff;

	return edge;
}


struct constraintRow
{
	vector<int> constraintNo;
};


void solve(int d, int N){

	vector<constraintRow> nodeConstraint;
	N = 2*N + 1;

	cout << "vertex_edge " << vertex_edge(5, 2, 3) << endl;

	int nodeCount = N;
	int maxN;
	maxN = (N * (N-1) * (N-2))/6 + N + (N * (N-1))/2;

	int constraint_number = (N * (N-1))/2;
	nodeConstraint.resize(constraint_number);

	std::chrono::steady_clock::time_point start, end;
	IloEnv env;
	IloModel model(env);
	IloBoolVarArray x = CreateBoolVarArray(env, (N * (N-1))/2, "x");

	IloObjective objective = IloMaximize(env);
	model.add(objective);

	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
		{
			objective.setLinearCoef(x[vertex_edge(N, i, j)], 1);
		}

	cout << "Constraint number initial " << constraint_number << endl;

	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			for (int k = j + 1; k < N; ++k)
			{
				model.add(x[vertex_edge(N, i, j)] + x[vertex_edge(N, i, k)] + x[vertex_edge(N, j, k)] <= 2);
				nodeConstraint[vertex_edge(N, i, j)].constraintNo.emplace_back(constraint_number);
				nodeConstraint[vertex_edge(N, i, k)].constraintNo.emplace_back(constraint_number);
				nodeConstraint[vertex_edge(N, j, k)].constraintNo.emplace_back(constraint_number);

				constraint_number++;
			}

	cout << "Constraint number after first set " << constraint_number << endl;

	IloRangeArray degreeConstraints = CreateRangeArray(env, N, "degree", d, d);

	for (int i = 0; i < N; ++i)
	{

		for (int j = 0; j < N; ++j)
			if (i < j)
			{
				degreeConstraints[i].setLinearCoef(x[vertex_edge(N, i, j)], 1);
				nodeConstraint[vertex_edge(N, i, j)].constraintNo.emplace_back(constraint_number);

			}
			else if (i > j)
			{
				degreeConstraints[i].setLinearCoef(x[vertex_edge(N, j, i)], 1);
				nodeConstraint[vertex_edge(N, j, i)].constraintNo.emplace_back(constraint_number);
			}
		constraint_number++;

	}

	cout << "Constraint number end " << constraint_number << endl;

	model.add(degreeConstraints);

	double maxUpperBound = (N * d) / 2;

	int not_maxed = 0;
	int currentBound = 0;
	int limitingBound = floor((N * d) / 2 - 0.5 * not_maxed);

	IloCplex cplex(model);
	int nThreads = 10;
	cplex.setParam(IloCplex::Threads, nThreads);
	//cplex.setParam(IloCplex::Param::TimeLimit, 120);

	//GlobalCallback CandidateCallback(currentBound);
	//BranchingCallback Callback(x, nodeConstraint, m, nodeCount, nThreads);
	//cplex.use(&Callback, IloCplex::Callback::Context::Id::Branching);
	//cplex.use(&CandidateCallback, IloCplex::Callback::Context::Id::Relaxation);

	not_maxed = 1;
	int timer = 0;
	int total_node_count = 0;
	int node_count = 0;
	int UB = 0;

	for (int j = 0; j < not_maxed; ++j) {
		degreeConstraints[j].setLB(0);
		degreeConstraints[j].setUB(d - 1);
	}

	start = chrono::steady_clock::now();

	bool Success = cplex.solve();
	if (not Success) {
		maxUpperBound = floor((N * d) / 2 - 0.5 * not_maxed);
		node_count = 0;
	}else{
		node_count = cplex.getNnodes();
	}

	if ((Success) && (cplex.getObjValue() > currentBound)){
		currentBound = cplex.getObjValue();
	}


	if ((Success) && (cplex.getObjValue() > currentBound)) {
		currentBound = cplex.getObjValue();
	}

	if (Success){
		cout << "For " << not_maxed << " edges; upper bound is " << cplex.getBestObjValue() << " lower bound is " << cplex.getObjValue() << endl;
	}

	cout << "For " << not_maxed << " edges; current bound is " << currentBound << endl;


	using namespace boost;
	{
		for (int vertex_no=0; vertex_no < N; ++vertex_no){
			typedef adjacency_list <vecS, vecS, undirectedS> Graph;
			Graph G;

			std::vector<graph_traits<Graph>::vertex_descriptor > mate(N);

			for (int i = 0; i < N; ++i)
				for (int j = i + 1; j < N; ++j) {
					IloNum xVal;
					xVal = cplex.getValue(x[vertex_edge(N, i, j)]);
					if (xVal > 0.5) {
						add_edge(i, j, G);
					}
				}

			remove_vertex(vertex_no, G);

			bool success = checked_edmonds_maximum_cardinality_matching(G, &mate[0]);

			int match_size = matching_size(G, &mate[0]);
			cout << "Removed vertex " << vertex_no << ", match size" << match_size << endl;
		}

	}

	end = chrono::steady_clock::now();
	int secs = chrono::duration_cast<chrono::seconds>(end - start).count();
	timer = timer + secs;


	cplex.getObjective();
	cplex.exportModel ("lpex1.lp");

	cout << "total time " << timer << endl;
	cout << "node count " << total_node_count << endl;
	cout << "best solution " << currentBound << endl;


	ofstream MyFile("results_combined_x_10.txt", std::ios_base::app);

	int UB1 = floor(currentBound);
	int UB2 = floor(cplex.getBestObjValue() + 0.01);

	int finalUB = max(UB1, UB2);
	finalUB = max(UB, finalUB);
	cout << "UB is " << UB2 << " but should be " << cplex.getBestObjValue() << endl;
	MyFile << d << "," << N << "," << currentBound << "," << finalUB << "," << timer << "," << total_node_count << endl;
	MyFile.close();

}

int main() {

	ofstream MyFile("results_combined_x_10.txt");
	MyFile << "d,N,LB,UB,time,nodeCount \n";
	MyFile.close();

	vector<vector<int>> params {
	/*	{7,	8},
		{7,	9},
		{7, 10},
		{8,	9},
		{8, 10},
		{9, 10},
		{9, 11},
		{9, 12},
		 {10, 11},
		{10, 12},
		{11, 12},
		{11, 13},
		{11, 14},
		{11, 15},
		{12, 13},
		{12, 14},
		{12, 15}, */
		//{13, 14},
		{11, 13},
		//{12, 14},
		//{13, 14},
		//{13, 15}
	};
	int d = 7;
	int N = 15;

	cout << params.size() << endl;

	for (int i=0; i<params.size(); i++){
		solve(params[i][0], params[i][1]);
	}

	return 0;
};
