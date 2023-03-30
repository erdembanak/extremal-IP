#include <iostream>
#include <chrono>
#include <Common.h>
#include <vector>
#include <mutex>
#include <typeinfo>

#define MAXN 8000
#define WORDSIZE 64

#define REMOVEONEARC0(g,v,w,m) DELELEMENT0(GRAPHROW0(g,v,m),w)
#define REMOVEONEEDGE(g,v,w,m) { REMOVEONEARC0(g,v,w,m); REMOVEONEARC0(g,w,v,m); }

mutex theMutex;

extern "C" {
	#include "nauty.h"
	#include "nausparse.h"
}

int vertex_edge(int N, int i, int j){
	int upper_lim = N - 1;
	int lower_lim = N - i;
	int diff = j - i - 1;

	int edge = (i * (upper_lim + lower_lim)) / 2 + diff;

	/*
	if(i == 0){
		edge = j - 1;
	}else{
		edge = (i * (upper_lim + lower_lim)) / 2  + diff;
	}

	*/
	return edge;
}

vector<int> branch_indexes(int orbits[], int size_m, vector<int> emptyNodes){

	int n[size_m] = {};

	for(int i=0; i<size_m; ++i){
		n[orbits[i]]++;
	}

	multimap<int, int> counts;

	for (int i=0; i<size_m; ++i){
		if (n[i] > 0){
			counts.insert(pair<int, int>(n[i], i));
		}
	}

	int main_var;
	for(auto it = counts.rbegin(); it != counts.rend(); ++it)
	{
		if (binary_search(emptyNodes.begin(), emptyNodes.end(), it->second) == 0){
			main_var = it->second;
			break;
		}
	}

	vector<int> branch_on;

	for (int i=0; i<size_m; ++i){
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
			ADDONEEDGE(g,i,j,m);
		}
	}
};

struct BranchingWorker
{
	vector<constraintRow>& nodeConstraint;
	int& m;
	BranchingWorker(vector<constraintRow>& nodeConstraint, int& m);
	graph g[MAXN*MAXM];

	vector<int> Solve(int n, vector<int> arrayRemove);
};

BranchingWorker::BranchingWorker(vector<constraintRow>& nodeConstraintIn, int& mIn) : nodeConstraint(nodeConstraintIn), m(mIn)
{
};

vector<int> BranchingWorker::Solve(int n, const vector<int> arrayRemove){
	int lab[n],ptn[n],orbits[n];
	DEFAULTOPTIONS_GRAPH(options);

	statsblk stats;

	EMPTYGRAPH(g,m,n);
	graphFromConstraintMatrix(g, m, nodeConstraint);

	int selectedNode;
	for (int i=0; i<arrayRemove.size(); ++i){
		selectedNode = arrayRemove[i];
		for (int j=0; j<nodeConstraint[selectedNode].constraintNo.size(); ++j){
			REMOVEONEEDGE(g,selectedNode,nodeConstraint[selectedNode].constraintNo[j],m);
		}
	}

	densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);
	cout << "Finish solve" << endl;

	int dummy_array[n] = {};

	for(int i=0; i<n; ++i){
		dummy_array[orbits[i]]++;
	}

	multimap<int, int> counts;

	for (int i=0; i<n; ++i){
		if (dummy_array[i] > 0){
			counts.insert(pair<int, int>(dummy_array[i], i));
		}
	}

	int main_var;
	for(auto it = counts.rbegin(); it != counts.rend(); ++it)
	{
		if (binary_search(arrayRemove.begin(), arrayRemove.end(), it->second) == 0){
			if(it->second < n){
				main_var = it->second;
				break;
			}
		}
	}

	vector<int> branch_on;

	for (int i=0; i<n; ++i){
		if (orbits[i] == main_var){
			branch_on.emplace_back(i);
		}
	}

	return branch_on;
}

struct BranchingCallback : public IloCplex::Callback::Function
{
	vector<BranchingWorker*> workers;
	IloBoolVarArray x;
	vector<constraintRow> nodeConstraint;
	int m;
	int nThreads;

	BranchingCallback(IloBoolVarArray& x, vector<constraintRow>& nodeConstraint, int& m, int& nThreads);

	void invoke(const IloCplex::Callback::Context& context);
};

BranchingCallback::BranchingCallback(
	IloBoolVarArray& xIn,
	vector<constraintRow>& nodeConstraintIn,
	int& mIn, int& nThreadsIn) : x(xIn), nodeConstraint(nodeConstraintIn), m(mIn), nThreads(nThreadsIn)
{
	workers.resize(nThreads);
	for (int i = 0; i < nThreads; ++i)
	{
		workers[i] = new BranchingWorker(nodeConstraintIn, m);
	};
}


void BranchingCallback::invoke(const IloCplex::Callback::Context& context)
{
	cout << "BranchingCallback::invoke has been called" << endl;
	int n = nodeConstraint.size();

	int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);
	cout << "thread no " << threadNo << endl;
	BranchingWorker* pWorker = workers[threadNo];

	double rel_obj = context.getRelaxationObjective();

	IloEnv env = context.getEnv();
	double obj = context.getRelaxationObjective();

	IloNumArray lowerBound(env);
	context.getLocalLB(x, lowerBound);

	IloNumArray upperBound(env);
	context.getLocalUB(x, upperBound);

	vector<int> arrayRemove;

	for(int i=0; i<x.getSize(); ++i){
		if(upperBound[i] == 0){
			arrayRemove.emplace_back(i);
		}
	}

	for(int i=0; i<x.getSize(); ++i){
		if(lowerBound[i] == 1){
			arrayRemove.emplace_back(i);
		}
	}
	cout << "Start solve" << endl;

	vector<int> branch_on = pWorker->Solve(n, arrayRemove);

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
	cout << "create branch " << endl;

}


int main() {

	vector<constraintRow> nodeConstraint;
	cout << "vertex_edge " << vertex_edge(5, 2, 3) << endl;
	int d = 7;
	int N = 15;
	int maxN;
	maxN = (N * (N-1) * (N-2))/6 + N + (N * (N-1))/2;

	int n,m;
	n = maxN;

	m = ceil(n/WORDSIZE);
	m = SETWORDSNEEDED(n);

	graph g[MAXN*MAXM];
	EMPTYGRAPH(g,m,n);

	cout << g << endl;
	cout << *g << endl;

	/*
	/*
	graph g[MAXN*MAXM];
	set h[MAXN*MAXM];
	int lab[maxN],ptn[maxN],orbits[maxN];
	static DEFAULTOPTIONS_GRAPH(options);

	statsblk stats;


	EMPTYGRAPH(h,m,n);
	*/


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
				//ADDONEEDGE(g,vertex_edge(N, i, j),constraint_number,m);
				//ADDONEEDGE(g,vertex_edge(N, i, k),constraint_number,m);
				//ADDONEEDGE(g,vertex_edge(N, j, k),constraint_number,m);
				nodeConstraint[vertex_edge(N, i, j)].constraintNo.emplace_back(constraint_number);
				nodeConstraint[vertex_edge(N, i, k)].constraintNo.emplace_back(constraint_number);
				nodeConstraint[vertex_edge(N, j, k)].constraintNo.emplace_back(constraint_number);

				//REMOVEONEEDGE(g,vertex_edge(N, i, j),constraint_number,m);
				//REMOVEONEEDGE(g,vertex_edge(N, i, k),constraint_number,m);
				//REMOVEONEEDGE(g,vertex_edge(N, j, k),constraint_number,m);

				//cout << "Edge btw " << constraint_number << " and " << vertex_edge(N, i, j) << endl;
				//cout << "Edge btw " << constraint_number << " and " << vertex_edge(N, i, k) << endl;
				//cout << "Edge btw " << constraint_number << " and " << vertex_edge(N, j, k) << endl;
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
				//ADDONEEDGE(g, vertex_edge(N, i, j), constraint_number,m);
				nodeConstraint[vertex_edge(N, i, j)].constraintNo.emplace_back(constraint_number);

				//cout << "Edge btw " << constraint_number << " and " << vertex_edge(N, i, j) << endl;
				//REMOVEONEEDGE(g, vertex_edge(N, i, j), constraint_number,m);
			}
			else if (i > j)
			{
				degreeConstraints[i].setLinearCoef(x[vertex_edge(N, j, i)], 1);
				//ADDONEEDGE(g,vertex_edge(N, j, i), constraint_number,m);
				nodeConstraint[vertex_edge(N, j, i)].constraintNo.emplace_back(constraint_number);

				//cout << "Edge btw " << constraint_number << " and " << vertex_edge(N, j, i) << endl;
				//REMOVEONEEDGE(g, vertex_edge(N, j, i), constraint_number,m);

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

/*
	graphFromConstraintMatrix(h, m, nodeConstraint);

	densenauty(h,lab,ptn,orbits,&options,&stats,m,n,NULL);

	end = chrono::steady_clock::now();
	int secs = chrono::duration_cast<chrono::milliseconds>(end - start).count();

	cout << "Duration is " << secs << endl;
	vector<int> vec;

	for(int i=0; i<maxN; ++i){
		vec.emplace_back(orbits[i]);
	}
	sort( vec.begin(), vec.end() );
	vec.erase( unique( vec.begin(), vec.end() ), vec.end() );

	cout << "size is " << vec.size() << endl;
*/

	int nThreads = 10;
	cplex.setParam(IloCplex::Threads, nThreads);

	BranchingCallback Callback(x, nodeConstraint, m, nThreads);
	cplex.use(&Callback, IloCplex::Callback::Context::Id::Branching);
	bool Success = cplex.solve();
	if (Success)
	{
		cout << "Success " << cplex.getCplexStatus() << endl;
		cout << "LB = " << cplex.getObjValue() << endl;
		cout << "UB = " << cplex.getBestObjValue() << endl;
		cout << "node " << cplex.getNnodes() << endl;
	}

	return 0;
};
