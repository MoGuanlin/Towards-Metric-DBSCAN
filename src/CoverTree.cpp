#include "../header/CoverTree.h"
using namespace std;

bool CoverTree::Insert(Point *p, const vector<Node *> &Q_layer, int layer, vector<vector<double>> &distance_matrix)
{

    vector<Node *> Q = Q_layer;
    float d_p_Q_layer = MAX_DIST;
    float d_p_Q = MAX_DIST;
    Node *closest_node;
    float max_child_dist = pow(2, layer);
    for (auto node : Q_layer)
    {

        if (node->children.count(layer - 1) > 0)
        {
            for (auto child : node->children[layer - 1])
            {
                Q.push_back(child);

                double dist;
                if (distance_matrix[p->ID][child->point_ptr->ID] == -1)
                    dist = distance_matrix[p->ID][child->point_ptr->ID] = distance_matrix[child->point_ptr->ID][p->ID] = NormDist(*p, *(child->point_ptr));
                else
                    dist = distance_matrix[p->ID][child->point_ptr->ID];

                if (dist < d_p_Q)
                    d_p_Q = dist;
            }
        }
    }

    for (auto node : Q_layer)
    {
        double dist;
        if (distance_matrix[p->ID][node->point_ptr->ID] == -1)
            dist = distance_matrix[p->ID][node->point_ptr->ID] = distance_matrix[node->point_ptr->ID][p->ID] = NormDist(*p, *(node->point_ptr));
        else
            dist = distance_matrix[p->ID][node->point_ptr->ID];
        if (dist < d_p_Q_layer)
        {
            d_p_Q_layer = dist;
            closest_node = node;
        }
    }

    d_p_Q = min(d_p_Q, d_p_Q_layer);

    if (d_p_Q > max_child_dist)
        return false;

    if (d_p_Q < MIN_DIST)
    {
        cout << "重复的point!" << endl;
        return true;
    }

    vector<Node *> next_Q_layer;
    for (auto node : Q)
    {
        if (distance_matrix[p->ID][node->point_ptr->ID] <= max_child_dist)
            next_Q_layer.push_back(node);
    }
    bool flag = this->Insert(p, next_Q_layer, layer - 1, distance_matrix);
    if (flag == false && d_p_Q_layer <= max_child_dist)
    {
        Node *p_node = new Node(layer - 1, p);
        p_node->parent = closest_node;
        closest_node->children[layer - 1].push_back(p_node);
        return true;
    }
    else if (flag == true)
    {
        return true;
    }
    else
        return false;
}

void CoverTree::Construct(vector<Point *> &data, vector<vector<double>> &distance_matrix)
{

    this->root = new Node(this->max_layer, data[0]);
    int n = data.size();
    vector<Node *> Q_max_layer;
    Q_max_layer.push_back(this->root);
    for (int i = 1; i < n; i++)
    {
        this->Insert(data[i], Q_max_layer, max_layer, distance_matrix);
    }
}

CoverTree::CoverTree(int max_l, vector<Point *> &data, vector<vector<double>> &distance_matrix)
{
    // this->root = new Node();
    this->max_layer = max_l;
    this->diameter_upperbound = pow(2, max_layer);
    this->Construct(data, distance_matrix);
};