#include "Common.h"
#include <fstream>
#include <chrono>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/bimap/bimap.hpp>

struct Component
{
int ID;
int Match;
int Edge;

Component(int IDIn, int MatchIn, int EdgeIn)
{
ID = IDIn;
Match = MatchIn;
Edge = EdgeIn;
}
};

void solveModel(int MaxMatch)
{
std::chrono::steady_clock::time_point start, end;
int noThreads = 6;
IloEnv env;
IloModel model(env);
vector<Component> Components;

Components.emplace_back(0, 1, 11);
Components.emplace_back(1, 11, 121);
Components.emplace_back(2, 12, 134);
Components.emplace_back(3, 13, 146);

IloIntVarArray x = CreateIntVarArray(env, Components.size(), "x", 0, MaxMatch);

IloObjective objective = IloMaximize(env);

for (int i = 0; i < Components.size(); ++i) {
objective.setLinearCoef(x[i], Components[i].Edge);
}

model.add(objective);

IloExpr expr = IloExpr(env);
for (int i = 0; i < Components.size(); ++i) {
expr += x[i] * Components[i].Match;
}

model.add(expr <= MaxMatch);

IloCplex cplex(model);
cplex.exportModel("knapsack.lp");

bool success = cplex.solve();

cout << "LB = " << cplex.getObjValue() << endl;

for (int i = 0; i < Components.size(); ++i)
{
cout << "From " << i << " used " << cplex.getValue(x[i]) << endl;
}

ofstream MyFile("knapsack_results_11.txt", std::ios_base::app);

MyFile << MaxMatch << "," << cplex.getObjValue() << "," << cplex.getValue(x[0]) << "," << cplex.getValue(x[1]) << "," << cplex.getValue(x[2]) << "," << cplex.getValue(x[3]) << endl;
MyFile.close();

}

int main()
{
ofstream MyFile("knapsack_results_11.txt");
MyFile << "maxMatch,obj, comp_1, comp_11, comp_12, comp_13 \n";
MyFile.close();

for (int i = 15; i < 500; ++i) {
solveModel(i);
}

return 0;
}

--
Ali Erdem Banak
Invent Analytics ⋰ Data Scientist
ARI-2 Teknokent İTÜ Ayazağa Kampüsü, A Blok Kat:3 No: 302 Maslak, İstanbul
Tel: +90 (212) 286 1025 | Mobile: +90 539 799 00 39
erdem.banak@inventanalytics.com | www.inventanalytics.com
