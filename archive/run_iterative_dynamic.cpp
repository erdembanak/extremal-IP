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

#define MAXN 0
#define WORDSIZE 64

#define REMOVEONEARC0(g,v,w,m) DELELEMENT0(GRAPHROW0(g,v,m),w)
#define REMOVEONEEDGE(g,v,w,m) { REMOVEONEARC0(g,v,w,m); REMOVEONEARC0(g,w,v,m); }

mutex(theMutex);

extern "C" {
	#include "nauty.h"
	#include "nausparse.h"
}

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

vector<int> branch_indexes(int orbits[], int vertex_count, int variable_count, vector<int> emptyNodes){

	int orbit_array[vertex_count] = {};

	for(int i=0; i<vertex_count; ++i){
		orbit_array[orbits[i]]++;
	}

	multimap<int, int> counts;

	for (int i=0; i<variable_count; ++i){
		if (orbit_array[i] > 0){
			counts.insert(pair<int, int>(orbit_array[i], i));
		}
	}

	int main_var;
	for(auto it = counts.rbegin();it != counts.rend(); ++it)
	{
		if (it->second < variable_count)
			if (!count(emptyNodes.begin(), emptyNodes.end(), it->second)){
				main_var = it->second;
				//cout << "main var is " << main_var << endl;
				break;
			}
	}

	vector<int> branch_on;

	for (int i=0; i<variable_count; ++i){
		if (orbits[i] == main_var){
			branch_on.emplace_back(i);
		}
	}

	return branch_on;
}

struct constraintRow
{
	vector<int> constraintNo;
};

void graphFromConstraintMatrix(graph g[MAXN*MAXM], int m, vector<constraintRow> cm){
	for (int i=0; i<cm.size(); ++i){
		for (int j=0; j<cm[i].constraintNo.size(); ++j){
			ADDONEEDGE(g,i,cm[i].constraintNo[j],m);
		}
	}
}

struct BranchingWorker
{
	vector<constraintRow>& nodeConstraint;
	int& m;
	BranchingWorker(vector<constraintRow>& nodeConstraint, int& m);
	graph g[MAXN*MAXM];

	vector<int> Solve(int vertex_count, vector<int> arrayRemove);
};

BranchingWorker::BranchingWorker(vector<constraintRow>& nodeConstraintIn, int& mIn) : nodeConstraint(nodeConstraintIn), m(mIn)
{
};

struct BranchingCallback : public IloCplex::Callback::Function
{
	vector<BranchingWorker*> workers;
	IloBoolVarArray x;
	vector<constraintRow> nodeConstraint;
	int m;
	int nodeCount;
	int nThreads;

	BranchingCallback(IloBoolVarArray& x, vector<constraintRow>& nodeConstraint, int& m, int& nodeCount, int& nThreads);

	void invoke(const IloCplex::Callback::Context& context);
};

BranchingCallback::BranchingCallback(
		IloBoolVarArray& xIn,
		vector<constraintRow>& nodeConstraintIn,
		int& mIn,
		int& nodeCountIn,
		int& nThreadsIn) : x(xIn), nodeConstraint(nodeConstraintIn), m(mIn), nodeCount(nodeCountIn), nThreads(nThreadsIn)
{
	workers.resize(nThreads);
	for (int i = 0; i < nThreads; ++i)
	{
		workers[i] = new BranchingWorker(nodeConstraintIn, m);
	};
}

vector<int> BranchingWorker::Solve(int vertex_count, const vector<int> arrayRemove){

	int variable_count = nodeConstraint.size();

	DYNALLSTAT(graph,g,g_sz);
	DYNALLSTAT(int,lab,lab_sz);
	DYNALLSTAT(int,ptn,ptn_sz);
	DYNALLSTAT(int,orbits,orbits_sz);

	DYNALLOC2(graph,g,g_sz,m,vertex_count,"malloc");
	DYNALLOC1(int,lab,lab_sz,vertex_count,"malloc");
	DYNALLOC1(int,ptn,ptn_sz,vertex_count,"malloc");
	DYNALLOC1(int,orbits,orbits_sz,vertex_count,"malloc");

	{
		std::unique_lock<std::mutex> lock(theMutex);

		DEFAULTOPTIONS_GRAPH(options);
		statsblk stats;

		EMPTYGRAPH(g,m,vertex_count);
		graphFromConstraintMatrix(g, m, nodeConstraint);

		int selectedNode;
		for (int i=0; i<arrayRemove.size(); ++i){
			selectedNode = arrayRemove[i];
			for (int j=0; j<nodeConstraint[selectedNode].constraintNo.size(); ++j){
				REMOVEONEEDGE(g,selectedNode,nodeConstraint[selectedNode].constraintNo[j],m);
			}
		}

		densenauty(g,lab,ptn,orbits,&options,&stats,m,vertex_count,NULL);
	}

	return branch_indexes(orbits, vertex_count, variable_count, arrayRemove);

}


void BranchingCallback::invoke(const IloCplex::Callback::Context& context)
{

	int depth = context.getIntInfo(IloCplex::Callback::Context::Info::NodeDepth);

	if (depth > 3){
		return;
	}

	int variable_count = nodeConstraint.size();
	int N = nodeCount;
	int vertex_count = (N * (N-1) * (N-2))/6 + N + (N * (N-1))/2;
	int m = SETWORDSNEEDED(vertex_count);

	int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);

	BranchingWorker* pWorker = workers[threadNo];

	//cout << "BranchingCallback::invoke has been called" << endl;
	IloCplex::CplexStatus status = context.getRelaxationStatus(0);
	double obj = context.getRelaxationObjective();

	IloEnv env = context.getEnv();

	IloNumArray lowerBound(env);
	context.getLocalLB(x, lowerBound);

	IloNumArray upperBound(env);
	context.getLocalUB(x, upperBound);

	int lab[vertex_count],ptn[vertex_count],orbits[vertex_count];

	vector<int> arrayRemove;

	//cout << "upper bound is 0 for: ";
	for(int i=0; i<x.getSize(); ++i){
		if(upperBound[i] == 0){
			arrayRemove.emplace_back(i);
			//cout << i << " ";
		}
	}

	//cout << endl << "lower bound is 1 for: ";
	for(int i=0; i<x.getSize(); ++i){
		if(lowerBound[i] == 1){
			arrayRemove.emplace_back(i);
			//cout << i << " ";
		}
	}

	//cout << endl << "# of fixed variables: " << arrayRemove.size() << endl;

	vector<int> branch_on = pWorker->Solve(vertex_count, arrayRemove);

	double rel_obj = context.getRelaxationObjective();

	IloNumVarArray branch_vars(env);
	IloNumArray zeros(env);
	IloArray<IloCplex::BranchDirection> IloDown(env);

	for (int i=0; i<branch_on.size(); ++i){
		branch_vars.add(x[branch_on[i]]);
		zeros.add(0);
		IloDown.add(IloCplex::BranchDown);
	}
	IloArray<IloCplex::BranchDirection> down;
	CPXLONG upChild, downChild;

	upChild = context.makeBranch(branch_vars, zeros, IloDown, rel_obj);
	downChild = context.makeBranch(x[branch_on[0]], 1, IloCplex::BranchUp, rel_obj);

	(void)downChild;
	(void)upChild;

	branch_vars.end();
	IloDown.end();

	/*
	if (branch_on.size() > 1)
	{
		for (int i=0; i<branch_on.size(); ++i){
			branch_vars.add(x[branch_on[i]]);
			zeros.add(0);
			IloDown.add(IloCplex::BranchDown);
		}
		IloArray<IloCplex::BranchDirection> down;
		CPXLONG upChild, downChild;

		upChild = context.makeBranch(branch_vars, zeros, IloDown, rel_obj);
		downChild = context.makeBranch(x[branch_on[0]], 1, IloCplex::BranchUp, rel_obj);

		(void)downChild;
		(void)upChild;

		branch_vars.end();
		IloDown.end();

	};
	*/
}

void solve(int d, int N){

	vector<constraintRow> nodeConstraint;
	N = 2*N + 1;

	cout << "vertex_edge " << vertex_edge(5, 2, 3) << endl;

	int nodeCount = N;
	int maxN;
	maxN = (N * (N-1) * (N-2))/6 + N + (N * (N-1))/2;
	//graph g[MAXN*MAXM];

	int n,m;
	n = maxN;

	m = ceil(n/WORDSIZE);
	m = SETWORDSNEEDED(n);
	nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

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
	cplex.setParam(IloCplex::Param::TimeLimit, 120);

	GlobalCallback CandidateCallback(currentBound);
	BranchingCallback Callback(x, nodeConstraint, m, nodeCount, nThreads);
	cplex.use(&Callback, IloCplex::Callback::Context::Id::Branching);
	//cplex.use(&CandidateCallback, IloCplex::Callback::Context::Id::Relaxation);

	not_maxed = 0;
	int timer = 0;
	int total_node_count = 0;
	int node_count = 0;
	int UB = 0;
	while (currentBound < limitingBound) {
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

		end = chrono::steady_clock::now();
		int secs = chrono::duration_cast<chrono::seconds>(end - start).count();
		timer = timer + secs;

		if (cplex.getBestObjValue() > currentBound) {
			cout << "Failed to find a good solution one by one" << endl;
			UB = floor(limitingBound - 0.5);
			break;
		}


        cplex.setParam(IloCplex::Param::TimeLimit, 1800 - timer);

		total_node_count = total_node_count + node_count;
		cout << "node count for " << not_maxed << " is " << node_count << endl;
 		if (timer > 1800){
			cout << "Failed to find a good solution one by one due to time limit" << endl;
			UB = floor(limitingBound - 0.5);
			break;
		}

		not_maxed++;
		limitingBound = floor((N * d) / 2 - 0.5 * not_maxed);

		if (currentBound > limitingBound){
			UB = currentBound;
		}
	}

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
		{12, 13},
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
