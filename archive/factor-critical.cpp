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

#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/adjacency_list.hpp>

#define MAXN 8000
#define WORDSIZE 64

mutex(theMutex);

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

	cplex.setParam(IloCplex::Param::TimeLimit, 1800);

	start = chrono::steady_clock::now();

	int nThreads = 10;
	cplex.setParam(IloCplex::Threads, nThreads);

	bool Success = cplex.solve();
	if (Success)
	{
		cout << "Success " << cplex.getCplexStatus() << endl;
		cout << "LB = " << cplex.getObjValue() << endl;
		cout << "UB = " << cplex.getBestObjValue() << endl;
		cout << "node " << cplex.getNnodes() << endl;
	}

	end = chrono::steady_clock::now();
	int secs = chrono::duration_cast<chrono::seconds>(end - start).count();


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

	//ofstream MyFile("results_orbit.txt", std::ios_base::app);

	//MyFile << d << "," << N << "," << cplex.getObjValue() << "," << cplex.getBestObjValue() << "," << secs << "," << cplex.getNnodes() << endl;
	//MyFile.close();

}

int main() {

	//ofstream MyFile("results_orbit.txt");
	//MyFile << "d,N,LB,UB,time,nodeCount \n";
	//MyFile.close();

	vector<vector<int>> params {
		//{7,	9},
		//{9, 12},
		//{11, 15},
		{13, 17}
		/*
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
		{12, 15},
		{13, 14},
		{13, 15},
		{13, 16},
		{13, 17}
		*/
	};
	int d = 7;
	int N = 15;

	cout << params.size() << endl;

	for (int i=0; i<params.size(); i++){
		solve(params[i][0], params[i][1]);
	}

	return 0;
};
