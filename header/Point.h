#pragma once
#define MAX_DIST 0x3f3f3f3f
#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#include <iostream>
#include <vector>
#include <map>
#include <string.h>

using namespace std;
class Center;
class Point
{
public:
    int ID; // ID in dataset
    int clusterID;
    bool is_core;   // init false
    Center *center; // kcenter ball center
    vector<double> coords;
    string text;

    Point();
    Point(int d);
    Point(const Point &p);
    void set_coords(int idx, double value);
    double get_coords(int idx);
};

class Node
{
public:
    Point *point_ptr;
    Node *parent;
    // 对于每个node，要精确记录它在每一层的children，即每一个children要附带记录一个层数。以后插入新的点时，要根据层数检索它在那一层的children
    // 第i层的点包含layer=i的点和layer>i的所有点。第i层的children仅代表一个Node在第i层新增的children，不包含Node自身
    // 第i层的Q包含了：第i+1层的Q(i+1),和Q(i+1)的每一个点在第i层的children
    // vector<Node *> children;
    std::map<int, vector<Node *>> children;

    int layer;
    float min_separate_dist;
    Node();
    Node(int l, Point *point);
    Node(const Node &n);
};

class Center : public Point
{
public:
    int weight;
    bool queued;
    double r;
    vector<Point *> members;
    vector<Node *> tree_members;
    vector<Center *> Ap;
    Node *node_ptr;

    Center(const Point &p);
    Center(Node *n);
};

double NormDist(const Point &p1, const Point &p2);
double LevenshteinDist(const Point &p1, const Point &p2);