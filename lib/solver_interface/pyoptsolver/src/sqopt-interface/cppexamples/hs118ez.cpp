#include <stdio.h>
#include <string.h>
#include <iostream>

#include "sqoptEZ.h"

using namespace std;

int main(int argc, char **argv) {
  int n     = 2;
  int m     = 1;
  SqoptEZ prob(n, m);
  double Hval[2] = {1.0, 1.0};
  int Hrow[2] = {0, 1};
  int Hcol[3] = {0, 1, 2};
  prob.setMatrixCSC('H', Hval, Hrow, Hcol);
  prob.setLinearCostToZero();
  double Aval[2] = {1.0, 1.0};
  int Arow[2] = {0, 0};
  int Acol[3] = {0, 1, 2};
  prob.setMatrixCSC('A', Aval, Arow, Acol);
  double lb[3] = {-4, -4, 1};
  double ub[3] = {4, 4, 1};
  prob.setVector('l', lb);
  prob.setVector('u', ub);
  int info = prob.solve(0);
  std::cout << "cost is " << prob.getObjective();
  return 1;
}
