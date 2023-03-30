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
#include <iterator>
#include <numeric>

using namespace std;

#define MAXN 8000
#define WORDSIZE 64

#define REMOVEONEARC0(g,v,w,m) DELELEMENT0(GRAPHROW0(g,v,m),w)
#define REMOVEONEEDGE(g,v,w,m) { REMOVEONEARC0(g,v,w,m); REMOVEONEARC0(g,w,v,m); }

mutex(theMutex);

extern "C" {
	#include "nauty.h"
	#include "nausparse.h"
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

struct RHSrow
{
	int vertex_no;
	int no_variables;
	int RHS;
	int sign;
};

struct constraintInfo
{
	int RHS;
	int no_variables;
};

void graphFromConstraintMatrix(graph g[MAXN*MAXM], int m, vector<constraintRow> cm, int vertex_count){
	for (int i=0; i<cm.size(); ++i){
		for (int j=0; j<cm[i].constraintNo.size(); ++j){
			ADDONEEDGE(g,i,cm[i].constraintNo[j],m);
		}
	}
/*
	for (int i=136; i < 816; ++i){
		ADDONEEDGE(g,i,vertex_count+100,m);
	}

	for (int i=816; i < 833; ++i){
		ADDONEEDGE(g,i,vertex_count+99,m);
	}
*/
}

void update_value(map<int, vector<int>>& d, int key, int new_val) {
    map<int, vector<int>>::iterator it = d.find(key);
    if (it != d.end()) {
        it->second.push_back(new_val);
    }
    else {
         vector<int> v;
         v.push_back(new_val);
         d.insert(make_pair(key, v));
    }
}

struct BranchingWorker
{
	vector<constraintRow>& nodeConstraint;
	vector<RHSrow>& RHSRows;

	int& m;
	BranchingWorker(vector<constraintRow>& nodeConstraint, vector<RHSrow>& RHSRows, int& m);
	graph g[MAXN*MAXM];

	vector<int> Solve(int vertex_count, vector<RHSrow> RHSRows,vector<int> arrayRemove0, vector<int> arrayRemove1);
};

BranchingWorker::BranchingWorker(vector<constraintRow>& nodeConstraintIn, vector<RHSrow>& RHSRowsIn, int& mIn) : nodeConstraint(nodeConstraintIn), RHSRows(RHSRowsIn), m(mIn)
{
};

struct BranchingCallback : public IloCplex::Callback::Function
{
	vector<BranchingWorker*> workers;
	IloBoolVarArray x;
	vector<constraintRow> nodeConstraint;
	vector<RHSrow> RHSRows;
	int m;
	int nodeCount;
	int nThreads;

	BranchingCallback(IloBoolVarArray& x, vector<constraintRow>& nodeConstraint, vector<RHSrow>& RHSRows, int& m, int& nodeCount, int& nThreads);

	void invoke(const IloCplex::Callback::Context& context);
};

BranchingCallback::BranchingCallback(
		IloBoolVarArray& xIn,
		vector<constraintRow>& nodeConstraintIn,
		vector<RHSrow>& RHSRowsIn,
		int& mIn,
		int& nodeCountIn,
		int& nThreadsIn) : x(xIn), nodeConstraint(nodeConstraintIn), RHSRows(RHSRowsIn), m(mIn), nodeCount(nodeCountIn), nThreads(nThreadsIn)
{
	workers.resize(nThreads);
	for (int i = 0; i < nThreads; ++i)
	{
		workers[i] = new BranchingWorker(nodeConstraintIn, RHSRowsIn, m);
	};
}

vector<int> BranchingWorker::Solve(int vertex_count, vector<RHSrow> RHSRows, const vector<int> arrayRemove0, const vector<int> arrayRemove1){

	int lab[vertex_count],ptn[vertex_count],orbits[vertex_count];
	int variable_count = nodeConstraint.size();

	vector<RHSrow> local_RHSRows;
	for(int i=0; i<RHSRows.size(); ++i){
		local_RHSRows.emplace_back(RHSRows[i]);
	}

	{
		std::unique_lock<std::mutex> lock(theMutex);

		graph g[MAXN*MAXM];

		DEFAULTOPTIONS_GRAPH(options);
		options.defaultptn=FALSE;
		memset(ptn, 1, sizeof(ptn));

		for(int i=0; i<vertex_count; ++i){
			lab[i] = i;
		}

		statsblk stats;

		EMPTYGRAPH(g,m,vertex_count);

		graphFromConstraintMatrix(g, m, nodeConstraint, vertex_count);

		int selectedNode;
		int constraint_start = local_RHSRows[0].vertex_no;
		ptn[constraint_start-1] = 0;

		for (int i=0; i<arrayRemove0.size(); ++i){
			selectedNode = arrayRemove0[i];
			for (int j=0; j<nodeConstraint[selectedNode].constraintNo.size(); ++j){
				REMOVEONEEDGE(g,selectedNode,nodeConstraint[selectedNode].constraintNo[j],m);
				local_RHSRows[nodeConstraint[selectedNode].constraintNo[j] - constraint_start].no_variables--;
				//local_RHSRows[0].no_variables--;

			}
		}

		for (int i=0; i<arrayRemove1.size(); ++i){
			selectedNode = arrayRemove1[i];
			for (int j=0; j<nodeConstraint[selectedNode].constraintNo.size(); ++j){
				REMOVEONEEDGE(g,selectedNode,nodeConstraint[selectedNode].constraintNo[j],m);
				local_RHSRows[nodeConstraint[selectedNode].constraintNo[j] - constraint_start].no_variables--;
				local_RHSRows[nodeConstraint[selectedNode].constraintNo[j] - constraint_start].RHS--;
			}
		}

		map<int, vector<int>> RHS_var;

		for (int i=0; i < local_RHSRows.size(); ++i){
			int vertex_no = local_RHSRows[i].vertex_no;
			int RHS = local_RHSRows[i].RHS;
			update_value(RHS_var, RHS, vertex_no);
		}

		vector<int> zero_out;
		int fill_no = constraint_start;

		for (map<int, vector<int>>::iterator it = RHS_var.begin(); it != RHS_var.end(); it++){
			vector<int> RHS_vertex_list = it->second;
			for(int i=0; i<RHS_vertex_list.size(); ++i){
				if(i==RHS_vertex_list.size() - 1){
					zero_out.emplace_back(fill_no);
				}
				lab[fill_no] = RHS_vertex_list[i];
				fill_no++;
			}
		}

		for(int i=0; i<zero_out.size();++i){
			ptn[zero_out[i]] = 0;
		}


		for(int i=0; i<vertex_count; ++i){
			if(ptn[i] == 0){
			}
		}

		/*


		vector<int> changedRHS;
		vector<int> changedVar;

		int change_RHS = -1;

		for (map<int, vector<int>>::iterator it = RHS_var.begin(); it != RHS_var.end(); it++)
		{
			vector<int> check = it->second;
			sort(check.begin(), check.end());
			int uniqCnt = unique(check.begin(), check.end()) - check.begin();
			sort(check.begin(), check.end());


			vector<int>::iterator it2 = unique(check.begin(), check.end());

			check.resize(distance(check.begin(),it2));

			vector<int> check0;
			for(int i=0; i<check.size(); ++i){
				if (check[i] > 0){
					check0.emplace_back(check[i]);
				}
			}

			if (check0.size() > 1){
				for (int j=0; j < check.size(); ++j){
					changedRHS.emplace_back(check0[j]);
				}
			}
			changedVar.emplace_back(it->first);
		}
				int count = 0;


		vector<int> checkConstraints;

		for (int i=0; i < local_RHSRows.size(); ++i){
			int RHS_bind = local_RHSRows[i].RHS;
			int constraint_to_bind = local_RHSRows[i].vertex_no;
			int no_var = local_RHSRows[i].no_variables;

			//bool found = (find(changedRHS.begin(), changedRHS.end(), RHS_bind) != changedRHS.end());
			bool found = (find(changedVar.begin(), changedVar.end(), no_var) != changedVar.end());

			checkConstraints.emplace_back(constraint_to_bind);

			//ADDONEEDGE(g,constraint_to_bind, vertex_count + 10 - 5, m);
			//REMOVEONEEDGE(g,constraint_to_bind, vertex_count + 10 - 5, m);

			//ADDONEEDGE(g,constraint_to_bind, vertex_count, m);
			//REMOVEONEEDGE(g,constraint_to_bind, vertex_count, m);

			//cout << "bind " << vertex_count << endl;
			count++;

			if (found){
				//ADDONEEDGE(g,constraint_to_bind, vertex_count + RHS_bind, m);
				//ptn[constraint_to_bind] = RHS_bind;
				//count++;
				//cout << "For RHS " << RHS_bind <<  " changed" << endl;
			}
		}
		// order constraints in RHS values
		// then 1-1-1-1
		int uniqCnt = unique(changedVar.begin(), changedVar.end()) - changedVar.begin();
*/

		densenauty(g,lab,ptn,orbits,&options,&stats,m,vertex_count,NULL);
	}

	std::vector<int> arrayRemove = arrayRemove0;
	arrayRemove.insert(arrayRemove.end(), arrayRemove1.begin(), arrayRemove1.end());


	return branch_indexes(orbits, vertex_count, variable_count, arrayRemove);

}

void BranchingCallback::invoke(const IloCplex::Callback::Context& context)
{
	int depth = context.getIntInfo(IloCplex::Callback::Context::Info::NodeDepth);

	if (depth > 20){
		return;
	}

	//int variable_count = nodeConstraint.size();
	int N = nodeCount;
	int vertex_count = (N * (N-1) * (N-2))/6 + N + (N * (N-1))/2;
	//int m = SETWORDSNEEDED(vertex_count);

	int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);

	BranchingWorker* pWorker = workers[threadNo];

	//cout << "BranchingCallback::invoke has been called" << endl;
	//IloCplex::CplexStatus status = context.getRelaxationStatus(0);
	//double obj = context.getRelaxationObjective();

	IloEnv env = context.getEnv();

	IloNumArray lowerBound(env);
	context.getLocalLB(x, lowerBound);

	IloNumArray upperBound(env);
	context.getLocalUB(x, upperBound);

	//int lab[vertex_count],ptn[vertex_count],orbits[vertex_count];

	vector<int> arrayRemove0;
	vector<int> arrayRemove1;

	//cout << "upper bound is 0 for: ";
	for(int i=0; i<x.getSize(); ++i){
		if(upperBound[i] == 0){
			arrayRemove0.emplace_back(i);
			//cout << i << " ";
		}
	}

	//cout << endl << "lower bound is 1 for: ";
	for(int i=0; i<x.getSize(); ++i){
		if(lowerBound[i] == 1){
			arrayRemove1.emplace_back(i);
			//cout << i << " ";
		}
	}

	//cout << endl << "# of fixed variables: " << arrayRemove.size() << endl;

	vector<int> branch_on = pWorker->Solve(vertex_count, RHSRows, arrayRemove0, arrayRemove1);

	double rel_obj = context.getRelaxationObjective();

	IloNumVarArray branch_vars(env);
	IloNumArray zeros(env);
	IloArray<IloCplex::BranchDirection> IloDown(env);

	//cout << "size of orbit " << branch_on.size() << endl;
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
}

IloBoolVarArray x;
vector<constraintRow> nodeConstraint;
int m;
int nodeCount;

int main() {

	vector<constraintRow> nodeConstraint;
	vector<RHSrow> RHSRows;

	int d = 8;
	int N = 19;
	int nodeCount = N;
	int maxN;
	maxN = (N * (N-1) * (N-2))/6 + N + (N * (N-1))/2;
	//graph g[MAXN*MAXM];

	//int lab[maxN],ptn[maxN],orbits[maxN];
	//static DEFAULTOPTIONS_GRAPH(options);

	statsblk stats;
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

	int vertex_no;
	int no_variables;
	int RHS;
	int sign;

	for (int i = 0; i < N; ++i)
		for (int j = i + 1; j < N; ++j)
			for (int k = j + 1; k < N; ++k)
			{
				model.add(x[vertex_edge(N, i, j)] + x[vertex_edge(N, i, k)] + x[vertex_edge(N, j, k)] <= 2);
				nodeConstraint[vertex_edge(N, i, j)].constraintNo.emplace_back(constraint_number);
				nodeConstraint[vertex_edge(N, i, k)].constraintNo.emplace_back(constraint_number);
				nodeConstraint[vertex_edge(N, j, k)].constraintNo.emplace_back(constraint_number);

				RHSrow theRHS;
				theRHS.vertex_no = constraint_number;
				theRHS.no_variables = 3;
				theRHS.RHS = 2;
				theRHS.sign = 0;

				RHSRows.emplace_back(theRHS);
				constraint_number++;
			}

	cout << "Constraint number after first set " << constraint_number << endl;
	int first_set = constraint_number;

	IloRangeArray degreeConstraints = CreateRangeArray(env, N, "degree", 0, d);

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

		RHSrow theRHS;
		theRHS.vertex_no = constraint_number;
		theRHS.no_variables = N-1;
		theRHS.RHS = d;
		theRHS.sign = 0;

		RHSRows.emplace_back(theRHS);

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

	BranchingCallback Callback(x, nodeConstraint, RHSRows, m, nodeCount, nThreads);
	cplex.use(&Callback, IloCplex::Callback::Context::Id::Branching);

	bool Success = cplex.solve();
	if (Success)
	{
		cout << "Success " << cplex.getCplexStatus() << endl;
		cout << "LB = " << cplex.getObjValue() << endl;
		cout << "UB = " << cplex.getBestObjValue() << endl;
		cout << "node " << cplex.getNnodes() << endl;
	}

	/*
	not_maxed = 0;
	int timer = 0;
	int total_node_count = 0;
	int node_count = 0;

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

		cout << "For " << not_maxed << " edges; current bound is " << currentBound << endl;

		not_maxed++;
		limitingBound = floor((N * d) / 2 - 0.5 * not_maxed);

		if (cplex.getBestObjValue() > currentBound) {
			cout << "Failed to find a good solution one by one" << endl;
			break;
		}

		end = chrono::steady_clock::now();
		int secs = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		timer = timer + secs;
		total_node_count = total_node_count + node_count;
		cout << "node count for " << not_maxed << " is " << node_count << endl;
 		if (timer > 600000){
			cout << "Failed to find a good solution one by one due to time limit" << endl;

			break;
		}

	}

	cplex.getObjective();
	cplex.exportModel ("lpex1.lp");

	cplex.setParam(IloCplex::Param::TimeLimit, 6000);

	start = chrono::steady_clock::now();

	cout << "total time " << timer << endl;
	cout << "node count " << total_node_count << endl;
	cout << "best solution " << currentBound << endl;
*/
	/*
	int nThreads = 10;
	cplex.setParam(IloCplex::Threads, nThreads);

	BranchingCallback Callback(x, nodeConstraint, m, nodeCount, nThreads);
	cplex.use(&Callback, IloCplex::Callback::Context::Id::Branching);

	bool Success = cplex.solve();

	if (Success)
	{
		cout << "Success " << cplex.getCplexStatus() << endl;
		cout << "LB = " << cplex.getObjValue() << endl;
		cout << "UB = " << cplex.getBestObjValue() << endl;
		cout << "node " << cplex.getNnodes() << endl;
	}
*/

	//BranchingCallback Callback(x, nodeConstraint, m, nodeCount, nThreads);
	//cplex.use(&Callback, IloCplex::Callback::Context::Id::Branching);

	return 0;
}

