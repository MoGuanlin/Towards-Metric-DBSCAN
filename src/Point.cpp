#include "../header/Point.h"
#include <cmath>
using namespace std;

Point::Point() : center(NULL), is_core(false), clusterID(-1), text("") {}
Point::Point(int d) : center(NULL), coords(d), is_core(false), clusterID(-1), text("") {}
Point::Point(const Point &p) : coords(p.coords), center(p.center), ID(p.ID), is_core(p.is_core), clusterID(p.clusterID), text(p.text) {}
void Point::set_coords(int idx, double value)
{
    coords[idx] = value;
}
double Point::get_coords(int idx)
{
    return coords[idx];
}

Center::Center(const Point &p) : Point(p), weight(1), r(0), queued(false) {}
Center::Center(Node *n) : Point(*(n->point_ptr)), node_ptr(n), weight(0), r(0), queued(false) {}

Node::Node() : point_ptr(NULL), parent(NULL) {}
Node::Node(int l, Point *point) : point_ptr(point), parent(NULL), layer(l), min_separate_dist(pow(2, l)) {}
Node::Node(const Node &n) : point_ptr(n.point_ptr), parent(n.parent), children(n.children), layer(n.layer), min_separate_dist(n.min_separate_dist) {}

double NormDist(const Point &p1, const Point &p2)
{
    double sum = 0;
    int dim = p1.coords.size();
    for (int i = 0; i < dim; i++)
    {
        double diff = p1.coords[i] - p2.coords[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}
double
LevenshteinDist(const Point &p1, const Point &p2)
{
    /* cout << "*********p1: " << p1.text << endl;
    cout << "*********p2: " << p2.text << endl;*/
    int row = p1.text.length();
    int col = p2.text.length();

    int mat[row][col];

    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            if (i == 0)
            {
                mat[i][j] = j;
            }
            else if (j == 0)
            {
                mat[i][j] = i;
            }
            else
            {
                int cost = (p1.text[i - 1] == p2.text[j - 1]) ? 0 : 1;
                mat[i][j] = MIN3(mat[i - 1][j] + 1,
                                 mat[i][j - 1] + 1,
                                 mat[i - 1][j - 1] + cost);
            }
        }
    }
    // cout << "edit distance: " << mat[row - 1][col - 1] << endl;
    return mat[row - 1][col - 1];
}