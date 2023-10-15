#pragma once
#include "IO.h"
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <chrono>
#define MIN_DIST 0.00001

using namespace std;

class EpsilonNet
{
public:
    vector<Node *> centers;
    int layer;
    float min_separate_dist;
};

class CoverTree
{
public:
    Node *root;
    int max_layer;
    float diameter_upperbound;
    CoverTree(int max_l, vector<Point *> &data, vector<vector<double>> &distance_matrix);
    void Construct(vector<Point *> &data, vector<vector<double>> &distance_matrix);
    bool Insert(Point *p, const vector<Node *> &Q_layer, int layer, vector<vector<double>> &distance_matrix);
};
