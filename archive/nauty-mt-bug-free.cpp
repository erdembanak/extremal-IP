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

	int lab[vertex_count],ptn[vertex_count],orbits[vertex_count];
	int variable_count = nodeConstraint.size();

	{
		std::unique_lock<std::mutex> lock(theMutex);

		graph g[MAXN*MAXM];

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

ILOBRANCHCALLBACK4(MyBranch, IloBoolVarArray, x, vector<constraintRow>, nodeConstraint, int, m,
int, nodeCount){
	int variable_count = nodeConstraint.size();
	int N = nodeCount;
	int vertex_count = (N * (N-1) * (N-2))/6 + N + (N * (N-1))/2;
	int m = SETWORDSNEEDED(vertex_count);

	//cout << "BranchingCallback::invoke has been called" << endl;

	IloEnv env = getEnv();

	IloNumArray lowerBound(env);
	IloNumArray upperBound(env);

	getLBs(lowerBound, x);
	getUBs(upperBound, x);

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

	{
		std::unique_lock<std::mutex> lock(theMutex);

		graph g[MAXN*MAXM];

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

	vector<int> branch_on = branch_indexes(orbits, vertex_count, variable_count, arrayRemove);

	IloArray<IloCplex::BranchDirection> directions(env);
	IloNumArray bounds(env);

	double rel_obj = getObjValue();

	IloNumVarArray branch_vars(env);
	IloNumArray zeros(env);
	IloArray<IloCplex::BranchDirection> IloDown(env);

	/*
	for (int i=0; i<branch_on.size(); ++i){
		branch_vars.add(x[branch_on[i]]);
		zeros.add(0);
		IloDown.add(IloCplex::BranchDown);
	}
	IloArray<IloCplex::BranchDirection> down;

	//cout << "size of orbit is " << branch_on.size() << " and first element is " << branch_on[0] << endl;

	makeBranch(x[branch_on[0]], 1, IloCplex::BranchUp, getObjValue());
	makeBranch(branch_vars, zeros, IloDown, getObjValue());

	*/
	if (branch_on.size() > 1)
	{
		for (int i=0; i<branch_on.size(); ++i){
			branch_vars.add(x[branch_on[i]]);
			zeros.add(0);
			IloDown.add(IloCplex::BranchDown);
		}
		IloArray<IloCplex::BranchDirection> down;

		makeBranch(x[branch_on[0]], 1, IloCplex::BranchUp, getObjValue());
		makeBranch(branch_vars, zeros, IloDown, getObjValue());
	};

}

ILONODECALLBACK0(MySelect) {
   IloInt remainingNodes = getNremainingNodes();
   IloInt bestnode = -1;
   IloInt maxdepth = -1;
   IloInt mindepth = 999999;

   IloNum maxiisum = 0.0;
   for (IloInt i = 0; i < remainingNodes; i++) {
      IloInt depth = getDepth(i);
      IloNum iisum = getInfeasibilitySum(i);

      if ( (depth >= maxdepth)                   &&
           (depth > maxdepth || iisum > maxiisum)  ) {
         bestnode = i;
         maxdepth = depth;
         maxiisum = iisum;
      }

      /*
      if ( (depth <= mindepth)                   &&
           (depth < mindepth || iisum > maxiisum)  ) {
         bestnode = i;
         mindepth = depth;
         maxiisum = iisum;
      }
       */
   }
   if ( bestnode >= 0 ) selectNode(bestnode);
}

int main() {

	vector<constraintRow> nodeConstraint;

	cout << "vertex_edge " << vertex_edge(5, 2, 3) << endl;
	int d = 8;
	int N = 19;
	int nodeCount = N;
	int maxN;
	maxN = (N * (N-1) * (N-2))/6 + N + (N * (N-1))/2;
	//graph g[MAXN*MAXM];

	int lab[maxN],ptn[maxN],orbits[maxN];
	static DEFAULTOPTIONS_GRAPH(options);

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
		constraint_number++;

	}

	cout << "Constraint number end " << constraint_number << endl;

	model.add(degreeConstraints);

	IloCplex cplex(model);
	cplex.getObjective();
	cplex.exportModel ("lpex1.lp");

	cout << "Number of cols " << cplex.getNcols() << endl;
	cout << "Number of rows " << cplex.getNrows() << endl;

	cplex.setParam(IloCplex::Param::TimeLimit, 6000);

	start = chrono::steady_clock::now();

	int nThreads = 10;
	cplex.setParam(IloCplex::Threads, nThreads);

	//BranchingCallback Callback(x, nodeConstraint, m, nodeCount, nThreads);
	//cplex.use(&Callback, IloCplex::Callback::Context::Id::Branching);

	bool Success = cplex.solve();
	if (Success)
	{
		cout << "Success " << cplex.getCplexStatus() << endl;
		cout << "LB = " << cplex.getObjValue() << endl;
		cout << "UB = " << cplex.getBestObjValue() << endl;
		cout << "node " << cplex.getNnodes() << endl;
	}

	return 0;
}
 
