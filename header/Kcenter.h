#pragma once
#include "IO.h"
#include <cstdlib>
/* #include <Eigen/Core>
#include <Eigen/Dense> */
#include <time.h>
#include <chrono>
#include <algorithm>

using namespace std;
// using namespace Eigen;
struct dist_struct
{
    int ID;
    double *dists_ptr;
};
bool cmp(dist_struct a, dist_struct b);
void Kcenter(vector<Point *> &data, vector<Center *> &centers, double r, vector<vector<double>> &distance_matrix);
void Randomized_Kcenter(vector<Point *> &data, vector<Center *> &centers, double r, vector<vector<double>> &distance_matrix, int init_num, int farthest_num, int select_num);
void Kcenter_edit(vector<Point *> &data, vector<Center *> &centers, double r, vector<vector<double>> &distance_matrix); // void Kcenter_Eigen(vector<Point *> &data, vector<Center *> &centers, double r);
