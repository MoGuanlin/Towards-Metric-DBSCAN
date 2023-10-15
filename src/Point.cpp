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
    int row = p1.text.length(); /* p1 的长度 */
    int col = p2.text.length(); /* p2 的长度 */

    int mat[row][col]; /* C99 - variable-length array */

    for (int i = 0; i < row; ++i)
    { /* 数组的行 */
        for (int j = 0; j < col; ++j)
        { /* 数组的列 */
            if (i == 0)
            {
                mat[i][j] = j; /* 初始化第1行为 [ 0 1 2 ... ] */
            }
            else if (j == 0)
            {
                mat[i][j] = i; /* 初始化第1列为 [ 0 1 2 ... ] */
            }
            else
            {
                int cost = (p1.text[i - 1] == p2.text[j - 1]) ? 0 : 1; /* 记录s1[i-1]与s2[j-1]是否相等 */
                mat[i][j] = MIN3(mat[i - 1][j] + 1,                    /* 取三者的最小值 */
                                 mat[i][j - 1] + 1,
                                 mat[i - 1][j - 1] + cost);
            }
        }
    }
    // cout << "edit distance: " << mat[row - 1][col - 1] << endl;
    return mat[row - 1][col - 1];
}