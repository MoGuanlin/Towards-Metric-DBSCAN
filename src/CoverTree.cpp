#include "../header/CoverTree.h"
using namespace std;

bool CoverTree::Insert(Point *p, const vector<Node *> &Q_layer, int layer, vector<vector<double>> &distance_matrix)
{
    // 初始化Q = Q_layer中的所有点，外加Q_layer中每一个node在layer-1层的所有孩子
    vector<Node *> Q = Q_layer;
    float d_p_Q_layer = MAX_DIST;
    float d_p_Q = MAX_DIST;
    Node *closest_node; // Q_layer中与p距离最小的点
    float max_child_dist = pow(2, layer);
    for (auto node : Q_layer)
    {
        // 遍历Q_layer中的nodes，对于一个node，若它在第layer层有第layer-1层的孩子，则把这个孩子加入Q
        if (node->children.count(layer - 1) > 0)
        {
            for (auto child : node->children[layer - 1])
            {
                Q.push_back(child);
                // 顺便计算一下所有新加入的Q_layer的孩子与p之间的最短距离
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

    // 计算d(p,Q_layer)和d(p,Q)，并记录下Q_layer中与p最近的那个Node.

    // 计算d(p,Q_layer)
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

    // d(p,Q)为当前的d(p,Q)与d(p,Q_layer)二者取最小值
    d_p_Q = min(d_p_Q, d_p_Q_layer);

    // 如果p对于第i层及其孩子，都满足separate性质，说明p与第i层距离较远，至少和第i层节点是兄弟关系，不可能作为第i层以及第i层以下的节点的孩子
    if (d_p_Q > max_child_dist)
        return false;

    // 处理无法插入的重复point
    if (d_p_Q < MIN_DIST)
    {
        cout << "重复的point!" << endl;
        return true;
    }

    // 原文中的else
    //(a).将Q中与p距离小于等于2^i的点挑出来得到Q(i-1)
    vector<Node *> next_Q_layer;
    for (auto node : Q)
    {
        if (distance_matrix[p->ID][node->point_ptr->ID] <= max_child_dist)
            next_Q_layer.push_back(node);
    }
    //(b).如果向下递归查找parent失败，并且p可以作为第i层某点的孩子(即满足第i层的cover性质),则将p作为第i层的closest node的孩子，插入tree中
    bool flag = this->Insert(p, next_Q_layer, layer - 1, distance_matrix);
    if (flag == false && d_p_Q_layer <= max_child_dist)
    { // 在layer-1层为p创建一个node，将p作为closest node 的第layer-1层的孩子
        Node *p_node = new Node(layer - 1, p);
        p_node->parent = closest_node;
        closest_node->children[layer - 1].push_back(p_node); // 可以这样写吗?????????????????????????、
        return true;
    }
    else if (flag == true)
    { // 如果flag为真，说明向下递归查找parent成功，p已经被插入tree中，本层递归直接退出，返回true
        return true;
    }
    else // 如果向下递归查找parent失败，而且p也不能作为第i层某点的孩子(即与第i层节点相距太远，不满足第i层的cover性质)，
        // 那么p还有可能可以作为第i层上面的某层节点的孩子，本层递归返回false，去到上层查找parent。
        return false;
}

void CoverTree::Construct(vector<Point *> &data, vector<vector<double>> &distance_matrix)
{ // 建树过程：初始化root为根节点，从data[1]开始遍历，为每一个data[i]调用Insert,将其插入树中

    this->root = new Node(this->max_layer, data[0]); // data[0]作为root
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
    // 初始化：设定max_layer,diameter_upperbound,调用construct
    // this->root = new Node();
    this->max_layer = max_l;
    this->diameter_upperbound = pow(2, max_layer);
    this->Construct(data, distance_matrix);
};