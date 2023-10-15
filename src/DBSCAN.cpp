#include "../header/DBSCAN.h"
using namespace std;
void range_search(vector<Point *> &data, int MinPts, double epsilon, vector<vector<double>> distance_matrix)
{
    int n = data.size();
    int neighbor_count[n] = {0};
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            double dist = NormDist(*(data[i]), *(data[j]));
            distance_matrix[i][j] = distance_matrix[j][i] = dist;
            if (dist <= epsilon)
            {
                neighbor_count[i]++;
                if (i != j)
                    neighbor_count[j]++;
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (neighbor_count[i] >= MinPts)
            data[i]->is_core = true;
    }
}
void range_search_edit(vector<Point *> &data, int MinPts, double epsilon, vector<vector<double>> distance_matrix)
{
    int n = data.size();
    int neighbor_count[n] = {0};
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            double dist = LevenshteinDist(*(data[i]), *(data[j]));
            distance_matrix[i][j] = distance_matrix[j][i] = dist;
            if (dist <= epsilon)
            {
                neighbor_count[i]++;
                if (i != j)
                    neighbor_count[j]++;
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (neighbor_count[i] >= MinPts)
            data[i]->is_core = true;
    }
}
int cluster_core_vanilla(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix)
{
    queue<Point *> Q;
    int n = data.size();
    bool queued[n] = {0};
    int clusterID = 0;
    for (auto i : data)
    {
        if (!queued[i->ID])
        {
            queued[i->ID] = true;
            Q.push(i);
            while (!Q.empty())
            {
                Point *head = Q.front();
                head->clusterID = clusterID;
                Q.pop();
                for (auto j : data)
                {
                    if (distance_matrix[head->ID][j->ID] <= epsilon && j->is_core == true && !queued[j->ID])
                    {
                        queued[j->ID] = true;
                        Q.push(j);
                    }
                }
            }
            clusterID++;
        }
    }

    // cluster border
    for (auto i : data)
    {
        if (i->is_core == false)
        {
            for (auto j : data)
            {
                if (j->is_core == true && distance_matrix[i->ID][j->ID] <= epsilon)
                {
                    i->clusterID = j->clusterID;
                    break;
                }
            }
        }
    }
    return clusterID;
}
void calculate_epsilon_net(CoverTree *T, int epsilon_net_layer, vector<vector<double>> &distance_matrix, vector<Center *> &centers)

{ // 从root开始，遍历整棵树
    queue<Node *> Q;
    Q.push(T->root);
    while (!(Q.empty()))
    {
        Node *n = Q.front();
        Q.pop();
        // 把n每一层的所有孩子入队列
        for (auto iter : (n->children))
        {
            for (auto iter2 : (iter.second))
                Q.push(iter2);
        }
        // 若n位于epsilon_net_layer的上方，则n为epsilon net中的一个center
        if (n->layer >= epsilon_net_layer)
        {
            Center *new_center = new Center(n);
            centers.push_back(new_center);
        }
    }
}
void calculate_Ap(vector<Center *> &centers, double epsilon, vector<vector<double>> &distance_matrix)
{

    int center_num = centers.size();

    for (int i = 0; i < center_num; i++)
    { // calculate Ap
        for (int j = i; j < center_num; j++)
        {
            double dist = distance_matrix[centers[i]->ID][centers[j]->ID];
            if (dist <= centers[i]->r + centers[j]->r + epsilon)
            {
                centers[i]->Ap.push_back(centers[j]);
                if (i != j)
                    centers[j]->Ap.push_back(centers[i]);
            }
        }
    }
}
void calculate_Ap_covertree(vector<Center *> &centers, vector<Point *> &data, double r_max, double epsilon, vector<vector<double>> &distance_matrix)
{

    int center_num = centers.size();

    for (int i = 0; i < center_num; i++)
    {
        for (int j = i; j < center_num; j++)
        {
            double dist = distance_matrix[centers[i]->node_ptr->point_ptr->ID][centers[j]->node_ptr->point_ptr->ID];
            if (dist == -1)
            {
                dist = distance_matrix[centers[i]->node_ptr->point_ptr->ID][centers[j]->node_ptr->point_ptr->ID] = distance_matrix[centers[j]->node_ptr->point_ptr->ID][centers[i]->node_ptr->point_ptr->ID] = NormDist(*(data[centers[j]->node_ptr->point_ptr->ID]), *(data[centers[i]->node_ptr->point_ptr->ID]));
            }

            if (dist <= 2 * r_max + epsilon)
            {
                centers[i]->Ap.push_back(centers[j]);
                if (i != j)
                    centers[j]->Ap.push_back(centers[i]);
            }
        }
    }
}
void mark_core(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    int n = data.size();
    int center_num = centers.size();
    Center *new_center;
    int dense_count = 0;
    int sparse_count = 0;

    for (int i = 0; i < center_num; i++)
    {
        if (centers[i]->weight >= MinPts)
        { // if center i have more than MinPts members, all of the members are core points
            dense_count++;
            centers[i]->is_core = true;
            for (auto iter = centers[i]->members.begin(); iter != centers[i]->members.end(); iter++)
            {
                (*iter)->is_core = true;
            }
            C.push_back(centers[i]);
        }
        else
        {
            sparse_count++;
            // count neighbor points
            for (auto iter = centers[i]->members.begin(); iter != centers[i]->members.end(); iter++)
            {
                Point *p = *iter;
                int neighbor_count = 0;
                for (auto iter2 = centers[i]->Ap.begin(); iter2 != centers[i]->Ap.end(); iter2++)
                {
                    Center *c = *iter2;
                    for (auto iter3 = c->members.begin(); iter3 != c->members.end(); iter3++)
                    {
                        if (distance_matrix[p->ID][(*iter3)->ID] == -1)
                        {
                            distance_matrix[p->ID][(*iter3)->ID] = distance_matrix[(*iter3)->ID][p->ID] = NormDist(*p, **iter3);
                        }

                        if (distance_matrix[p->ID][(*iter3)->ID] <= epsilon)
                        {
                            neighbor_count++;
                        }
                    }
                }
                if (neighbor_count >= MinPts)
                { // add p to S_*
                    p->is_core = true;
                    new_center = new Center(*p);
                    p->center = new_center;
                    new_center->members.push_back(p);
                    new_center->is_core = true;
                    C.push_back(new_center);
                }
            }
        }
    }
    cout << "dense count: " << dense_count << endl;
    cout << "sparse count: " << sparse_count << endl;
}
void mark_core_edit(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    int n = data.size();
    int center_num = centers.size();
    Center *new_center;
    int dense_count = 0;
    int sparse_count = 0;

    for (int i = 0; i < center_num; i++)
    {
        if (centers[i]->weight >= MinPts)
        { // if center i have more than MinPts members, all of the members are core points
            dense_count++;
            centers[i]->is_core = true;
            for (auto iter = centers[i]->members.begin(); iter != centers[i]->members.end(); iter++)
            {
                (*iter)->is_core = true;
            }
            C.push_back(centers[i]);
        }
        else
        {
            sparse_count++;
            // count neighbor points
            for (auto iter = centers[i]->members.begin(); iter != centers[i]->members.end(); iter++)
            {
                Point *p = *iter;
                int neighbor_count = 0;
                for (auto iter2 = centers[i]->Ap.begin(); iter2 != centers[i]->Ap.end(); iter2++)
                {
                    Center *c = *iter2;
                    for (auto iter3 = c->members.begin(); iter3 != c->members.end(); iter3++)
                    {
                        if (distance_matrix[p->ID][(*iter3)->ID] == -1)
                        {
                            distance_matrix[p->ID][(*iter3)->ID] = distance_matrix[(*iter3)->ID][p->ID] = LevenshteinDist(*p, **iter3);
                        }

                        if (distance_matrix[p->ID][(*iter3)->ID] <= epsilon)
                        {
                            neighbor_count++;
                        }
                    }
                }
                if (neighbor_count >= MinPts)
                { // add p to S_*
                    p->is_core = true;
                    new_center = new Center(*p);
                    p->center = new_center;
                    new_center->members.push_back(p);
                    new_center->is_core = true;
                    C.push_back(new_center);
                }
            }
        }
    }
    cout << "dense count: " << dense_count << endl;
    cout << "sparse count: " << sparse_count << endl;
}
void mark_core_covertree(vector<Center *> &centers, CoverTree *T, vector<vector<double>> &distance_matrix, int minpts, double epsilon, vector<Point *> &data, vector<Center *> &C)
{ // 遍历centers，每一个center遍历子树，将每一个node加入center node的tree_members中，且center node weight+1
    // 第二次遍历centers，若一个center的weight大于minpts，则将其所有tree_members nodes的points都标记为core
    // 若一个center的weight小于minpts，则遍历其tree_members，每个member遍历Ap,对于每个Ap center，遍历它的tree_members进行距离查询，距自己epsilon以内的点多于minpts个则标记自身为core
    for (auto center : centers)
    {
        // center和自己的point_ptr同步
        center->center = center;
        // 遍历子树
        queue<Node *> Q;
        Node *center_node = center->node_ptr;
        Q.push(center_node);
        double max_dist = 0;
        while (!(Q.empty()))
        {
            Node *n = Q.front();
            Q.pop();
            // 统计当前center的tree_members和weight和r
            center->tree_members.push_back(n);
            center->members.push_back(n->point_ptr);
            center->weight += 1;

            double dist = distance_matrix[center_node->point_ptr->ID][n->point_ptr->ID];
            if (dist == -1)
                dist = distance_matrix[center_node->point_ptr->ID][n->point_ptr->ID] = distance_matrix[n->point_ptr->ID][center_node->point_ptr->ID] = NormDist(*(n->point_ptr), *(center_node->point_ptr));
            if (dist > max_dist)
            {
                max_dist = dist;
            }
            for (auto iter : (n->children))
            {
                for (auto iter2 : (iter.second))
                    Q.push(iter2);
            }
        }
        center->r = max_dist;
    }

    for (auto center : centers)
    {
        if (center->weight >= minpts) // dense center
        {
            for (auto member : center->tree_members)
            {
                member->point_ptr->is_core = true;
                member->point_ptr->center = center;
            }
            C.push_back(center);
        }

        else // sparse center
        {
            for (auto member : center->tree_members)
            {
                member->point_ptr->center = center;
                int neighbor_count = 0;
                for (auto ap : center->Ap) // 遍历每个Ap的所有tree_members
                {
                    for (auto neighbor_candidate : ap->tree_members)
                    {
                        double dist = distance_matrix[member->point_ptr->ID][neighbor_candidate->point_ptr->ID];
                        if (dist == -1)
                        {
                            dist = distance_matrix[member->point_ptr->ID][neighbor_candidate->point_ptr->ID] = dist = distance_matrix[neighbor_candidate->point_ptr->ID][member->point_ptr->ID] = NormDist(*(neighbor_candidate->point_ptr), *(member->point_ptr));
                        }
                        if (dist <= epsilon)
                        {
                            neighbor_count++;
                        }
                    }
                }
                if (neighbor_count >= minpts)
                {
                    member->point_ptr->is_core = true;
                    Center *new_center = new Center(*(member->point_ptr));
                    member->point_ptr->center = new_center;
                    new_center->members.push_back(member->point_ptr);
                    new_center->is_core = true;
                    C.push_back(new_center);
                }
            }
        }
    }
}
Node *NN_search(Node *p, Node *root, vector<vector<double>> &distance_matrix, vector<Point *> &data, int terminate_layer)
{ // 迭代Qi集合，如果满足条件的Qi-1为空，则在Qi中找离p最近的点，作为p的最近邻返回
    int layer = root->layer;
    vector<Node *> Qi;
    Qi.push_back(root);

    while (true)
    {
        vector<Node *> Q = Qi;
        // 找Qi中每个点在layer-1层新增的孩子，加入Qi
        for (auto node : Qi)
        {
            if (node->children.count(layer - 1) > 0)
                for (auto child : node->children[layer - 1])
                    Q.push_back(child);
        }
        // 筛选Q,选出那些距离p比较近的作为next_Qi
        vector<Node *> next_Qi;
        for (auto node : Q)
        {
            double dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID];
            if (dist == -1)
            {
                dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID] = dist = distance_matrix[node->point_ptr->ID][p->point_ptr->ID] = NormDist(*(node->point_ptr), *(p->point_ptr));
            }
            if (dist <= pow(2, layer))
            {
                next_Qi.push_back(node);
            }
        }
        if (next_Qi.size() == 0) // 如果没有next_Qi,说明Qi已经是p在root子树中的最近邻居集合了，需要在Qi中找一个距离p最近的点返回
        {
            double min_dist = MAX_DIST;
            Node *NN;
            for (auto node : Qi)
            {
                double dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID];
                if (dist == -1)
                {
                    dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID] = dist = distance_matrix[node->point_ptr->ID][p->point_ptr->ID] = NormDist(*(node->point_ptr), *(p->point_ptr));
                }
                if (dist <= min_dist)
                {
                    min_dist = dist;
                    NN = node;
                }
            }
            return NN;
        }
        else if (layer <= terminate_layer) // 如果layer已经小于阈值，则直接在当前的p最近邻候选集中寻找最近邻
        {
            double min_dist = MAX_DIST;
            Node *NN;
            for (auto node : next_Qi)
            {
                double dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID];
                if (dist == -1)
                {
                    dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID] = dist = distance_matrix[node->point_ptr->ID][p->point_ptr->ID] = NormDist(*(node->point_ptr), *(p->point_ptr));
                }
                if (dist <= min_dist)
                {
                    min_dist = dist;
                    NN = node;
                }
            }
            return NN;
        }
        else // 如果尚未到达layer的最小阈值，且next_Qi中仍然有node，则继续迭代
        {
            Qi = next_Qi;
            layer--;
        }
    }
}
Node *NN_search2(Node *p, Center *root, vector<vector<double>> &distance_matrix, vector<Point *> &data, int terminate_layer)
{ // 遍历root的每个member

    double min_dist = MAX_DIST;
    Node *NN;
    for (auto node : root->tree_members)
    {
        double dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID];
        if (dist == -1)
        {
            dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID] = dist = distance_matrix[node->point_ptr->ID][p->point_ptr->ID] = NormDist(*(node->point_ptr), *(p->point_ptr));
        }
        if (dist <= min_dist)
        {
            min_dist = dist;
            NN = node;
        }
    }
    return NN;
}
void connect_centers(vector<Center *> &centers, vector<vector<double>> &distance_matrix, vector<Point *> &data, int terminate_layer, double r_max, int minpts, double epsilon)
{
    auto start = chrono::high_resolution_clock::now();
    int cluster_ID = 0;
    queue<Center *> Q;
    if (centers.size() == 0)
    {
        cout << "center size == 0!" << endl;
        return;
    }
    for (auto head_center : centers)
    {
        if (head_center->queued == true)
            continue;
        if (head_center->queued == false)
            Q.push(head_center);
        head_center->queued == true;
        while (Q.empty() == false)
        {
            Center *center = Q.front();
            Q.pop();

            //  对center.member中每一个member p，先设置p的clusterID，然后如果p是core point，则遍历center的Aps：
            //  判断这个Ap的类别：
            //  如果Ap是dense，则调用NN_search寻找p的最近邻，看最近邻与p的距离,若这个距离小于等于epsilon，则将Ap入队，设置queued
            //  如果Ap是sparse，则遍历Ap的core tree_members，查询p与这些core tree_members的最小距离，若这个距离小于等于epsilon，则将Ap入队，设置queued
            for (auto p : center->tree_members)
            {

                if (p->point_ptr->is_core == true)
                {
                    p->point_ptr->clusterID = cluster_ID; // 如果p不是core，则等待之后的cluster border来确定p的clusterID

                    for (auto Ap : center->Ap)
                    {
                        if (Ap->weight >= minpts)
                        {
                            // Node *NN = NN_search(p, Ap->node_ptr, distance_matrix, data, terminate_layer);
                            Node *NN = NN_search2(p, Ap, distance_matrix, data, terminate_layer);
                            if (distance_matrix[p->point_ptr->ID][NN->point_ptr->ID] <= epsilon)
                            {
                                if (Ap->queued == false)
                                    Q.push(Ap);
                                Ap->queued = true;
                                continue;
                            }
                        }
                        else
                        {
                            Node *NN = NULL;
                            double min_dist = MAX_DIST;
                            for (auto node : Ap->tree_members)
                            {
                                if (node->point_ptr->is_core == true)
                                {
                                    double dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID];
                                    if (dist == -1)
                                    {
                                        dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID] = dist = distance_matrix[node->point_ptr->ID][p->point_ptr->ID] = NormDist(*(node->point_ptr), *(p->point_ptr));
                                    }
                                    if (dist <= min_dist)
                                    {
                                        min_dist = dist;
                                        NN = node;
                                    }
                                }
                            }
                            if (distance_matrix[p->point_ptr->ID][NN->point_ptr->ID] <= epsilon)
                            {
                                if (Ap->queued == false)
                                    Q.push(Ap);
                                Ap->queued = true;
                                continue;
                            }
                        }
                    }
                }
            }
        }
        cluster_ID++;
    }
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

    // cluster borders
    // 尚未标记所有的非core point p：对于每个非core p，遍历Ap的所有core point，遇到距离p epsilon以内的core，就把p的clusterID设为这个core的

    start = chrono::high_resolution_clock::now();
    for (auto center : centers)
    {
        if (center->tree_members.size() < minpts)
        {
            for (auto p : center->tree_members)
            {
                if (p->point_ptr->is_core == false)
                {
                    for (auto Ap : center->Ap)
                    {
                        if (Ap->weight >= minpts)
                        {
                            Node *NN = NN_search(p, Ap->node_ptr, distance_matrix, data, terminate_layer);
                            if (distance_matrix[p->point_ptr->ID][NN->point_ptr->ID] <= epsilon)
                            {
                                p->point_ptr->clusterID = NN->point_ptr->clusterID;
                                continue;
                            }
                        }
                        else
                        {
                            Node *NN = NULL;
                            double min_dist = MAX_DIST;
                            for (auto node : Ap->tree_members)
                            {
                                if (node->point_ptr->is_core == true)
                                {
                                    double dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID];
                                    if (dist == -1)
                                    {
                                        dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID] = dist = distance_matrix[node->point_ptr->ID][p->point_ptr->ID] = NormDist(*(node->point_ptr), *(p->point_ptr));
                                    }
                                    if (dist <= min_dist)
                                    {
                                        min_dist = dist;
                                        NN = node;
                                    }
                                }
                            }
                            if (distance_matrix[p->point_ptr->ID][NN->point_ptr->ID] <= epsilon)
                            {
                                p->point_ptr->clusterID = NN->point_ptr->clusterID;
                                continue;
                            }
                        }
                    }
                }
            }
        }
    }
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
}
void connect_centers_bruteforce(vector<Center *> &centers, vector<vector<double>> &distance_matrix, vector<Point *> &data, int terminate_layer, double r_max, int minpts, double epsilon)
{
    // connect center太慢！
    auto start = chrono::high_resolution_clock::now();
    int cluster_ID = 0;
    queue<Center *> Q;
    if (centers.size() == 0)
    {
        cout << "center size == 0!" << endl;
        return;
    }
    for (auto head_center : centers)
    {
        if (head_center->queued == true)
            continue;
        if (head_center->queued == false)
            Q.push(head_center);
        head_center->queued == true;
        while (Q.empty() == false)
        {
            Center *center = Q.front();
            Q.pop();

            //  对center.member中每一个member p，先设置p的clusterID，然后如果p是core point，则遍历center的Aps：
            //  判断这个Ap的类别：
            //  如果Ap是dense，则调用NN_search寻找p的最近邻，看最近邻与p的距离,若这个距离小于等于epsilon，则将Ap入队，设置queued
            //  如果Ap是sparse，则遍历Ap的core tree_members，查询p与这些core tree_members的最小距离，若这个距离小于等于epsilon，则将Ap入队，设置queued
            for (auto p : center->tree_members)
            {

                if (p->point_ptr->is_core == true)
                {
                    p->point_ptr->clusterID = cluster_ID; // 如果p不是core，则等待之后的cluster border来确定p的clusterID

                    for (auto Ap : center->Ap)
                    {
                        if (Ap->weight >= minpts)
                        {
                            // Node *NN = NN_search(p, Ap->node_ptr, distance_matrix, data, terminate_layer);
                            Node *NN = NN_search2(p, Ap, distance_matrix, data, terminate_layer);
                            if (distance_matrix[p->point_ptr->ID][NN->point_ptr->ID] <= epsilon)
                            {
                                if (Ap->queued == false)
                                    Q.push(Ap);
                                Ap->queued = true;
                                continue;
                            }
                        }
                        else
                        {
                            Node *NN = NULL;
                            double min_dist = MAX_DIST;
                            for (auto node : Ap->tree_members)
                            {
                                if (node->point_ptr->is_core == true)
                                {
                                    double dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID];
                                    if (dist == -1)
                                    {
                                        dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID] = dist = distance_matrix[node->point_ptr->ID][p->point_ptr->ID] = NormDist(*(node->point_ptr), *(p->point_ptr));
                                    }
                                    if (dist <= min_dist)
                                    {
                                        min_dist = dist;
                                        NN = node;
                                    }
                                }
                            }
                            if (distance_matrix[p->point_ptr->ID][NN->point_ptr->ID] <= epsilon)
                            {
                                if (Ap->queued == false)
                                    Q.push(Ap);
                                Ap->queued = true;
                                continue;
                            }
                        }
                    }
                }
            }
        }
        cluster_ID++;
    }
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Connect center time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

    // cluster borders
    // 尚未标记所有的非core point p：对于每个非core p，遍历Ap的所有core point，遇到距离p epsilon以内的core，就把p的clusterID设为这个core的

    start = chrono::high_resolution_clock::now();
    for (auto center : centers)
    {
        if (center->tree_members.size() < minpts)
        {
            for (auto p : center->tree_members)
            {
                if (p->point_ptr->is_core == false)
                {
                    bool cluster_found = false;
                    for (auto Ap : center->Ap)
                    { // 遍历Ap的所有core tree_members,遇到近的则break
                        for (auto node : Ap->tree_members)
                        {
                            if (node->point_ptr->is_core == true)
                            {
                                double dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID];
                                if (dist == -1)
                                {
                                    dist = distance_matrix[p->point_ptr->ID][node->point_ptr->ID] = dist = distance_matrix[node->point_ptr->ID][p->point_ptr->ID] = NormDist(*(node->point_ptr), *(p->point_ptr));
                                }
                                if (dist <= epsilon)
                                {
                                    p->point_ptr->clusterID = node->point_ptr->clusterID;
                                    cluster_found = true;
                                    break;
                                }
                            }
                        }
                        if (cluster_found)
                            break;
                    }
                }
            }
        }
    }
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
}
double BCP(vector<Point *> &c1, vector<Point *> &c2, double epsilon, vector<vector<double>> &distance_matrix)
{
    double min_dist = MAX_DIST;
    for (auto iter1 : c1)
    {
        for (auto iter2 : c2)
        {
            double dist = distance_matrix[iter1->ID][iter2->ID];
            if (dist == -1)
            {
                dist = distance_matrix[iter1->ID][iter2->ID] = distance_matrix[iter2->ID][iter1->ID] = NormDist((*iter1), (*iter2));
            }
            if (dist < min_dist)
            {
                min_dist = dist;
            }
        }
    }
    return min_dist;
}
double BCP_edit(vector<Point *> &c1, vector<Point *> &c2, double epsilon, vector<vector<double>> &distance_matrix)
{
    double min_dist = MAX_DIST;
    for (auto iter1 : c1)
    {
        for (auto iter2 : c2)
        {
            double dist = distance_matrix[iter1->ID][iter2->ID];
            if (dist == -1)
            {
                dist = distance_matrix[iter1->ID][iter2->ID] = distance_matrix[iter2->ID][iter1->ID] = LevenshteinDist((*iter1), (*iter2));
            }
            if (dist < min_dist)
            {
                min_dist = dist;
            }
        }
    }
    return min_dist;
}
int cluster_core(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    queue<Center *> Q;
    int cluster_num = 0;
    for (auto iter = C.begin(); iter != C.end(); iter++)
    {
        if ((*iter)->queued == false)
        {
            // cout << "start find Cluster " << cluster_num << endl;
            Q.push(*iter);
            (*iter)->queued = true;
            while (!(Q.empty()))
            {
                // dequeue
                Center *head_center = Q.front();
                head_center->clusterID = cluster_num;
                for (auto iter : head_center->members)
                {
                    iter->clusterID = cluster_num;
                }
                Q.pop();
                // find reachable neighbor in S*, add them to Queue
                for (auto Ap_iter = head_center->center->Ap.begin(); Ap_iter != head_center->center->Ap.end(); Ap_iter++)
                {

                    if ((*Ap_iter)->members.size() >= MinPts)
                    {
                        if ((*Ap_iter)->queued == false)
                        {

                            if (BCP(head_center->members, (*Ap_iter)->members, epsilon, distance_matrix) <= epsilon)
                            {
                                Q.push(*Ap_iter);
                                (*Ap_iter)->queued = true;
                            }
                        }
                    }
                    else
                    {

                        for (auto ApMember_iter = (*Ap_iter)->members.begin(); ApMember_iter != (*Ap_iter)->members.end(); ApMember_iter++)
                        {
                            if ((*ApMember_iter)->is_core == true)
                            {

                                Center *core_center = (*ApMember_iter)->center;
                                if (core_center->queued == false)
                                {
                                    if (BCP(head_center->members, core_center->members, epsilon, distance_matrix) <= epsilon)
                                    {
                                        Q.push(core_center);
                                        core_center->queued = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            cluster_num++;
        }
    }
    cout << "Max cluster number is : " << cluster_num << endl;
    return cluster_num;
}
int cluster_core_approx(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    queue<Center *> Q;
    int cluster_num = 0;
    for (auto iter = C.begin(); iter != C.end(); iter++)
    {
        if ((*iter)->queued == false)
        {
            // cout << "start find Cluster " << cluster_num << endl;

            Center *test_c = (*iter);
            Q.push(*iter);
            (*iter)->queued = true;
            while (!(Q.empty()))
            {
                // dequeue
                Center *head_center = Q.front();
                head_center->clusterID = cluster_num;
                Q.pop();
                // find reachable neighbor in S*, add them to Queue
                /* int member_count = 0; */
                for (auto member_iter = head_center->members.begin(); member_iter != head_center->members.end(); member_iter++)
                {
                    /* member_count++; */
                    /*                     cout << "member_iter is the " << member_count << " th member." << endl;
                                        cout << "head center has " << head_center->members.size() << " members." << endl; */

                    (*member_iter)->clusterID = cluster_num;

                    for (auto Ap_iter = head_center->center->Ap.begin(); Ap_iter != head_center->center->Ap.end(); Ap_iter++)
                    {

                        if ((*Ap_iter)->members.size() >= MinPts)
                        {
                            if ((*Ap_iter)->queued == false)
                            {
                                double dist = distance_matrix[(*member_iter)->ID][(*Ap_iter)->ID];
                                if (dist != -1 && dist <= (*Ap_iter)->r + epsilon)
                                {
                                    Q.push(*Ap_iter);
                                    (*Ap_iter)->queued = true;
                                }
                            }
                        }
                        else
                        {

                            for (auto ApMember_iter = (*Ap_iter)->members.begin(); ApMember_iter != (*Ap_iter)->members.end(); ApMember_iter++)
                            {
                                if ((*ApMember_iter)->is_core == true)
                                {

                                    Center *core_center = (*ApMember_iter)->center;
                                    if (core_center->queued == false)
                                    {
                                        if (distance_matrix[core_center->ID][(*member_iter)->ID] <= epsilon)
                                        {
                                            Q.push(core_center);
                                            core_center->queued = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            cluster_num++;
        }
    }
    cout << "Max cluster number is : " << cluster_num << endl;
    return cluster_num;
}
int cluster_core_edit(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    queue<Center *> Q;
    int cluster_num = 0;
    for (auto iter = C.begin(); iter != C.end(); iter++)
    {
        if ((*iter)->queued == false)
        {
            // cout << "start find Cluster " << cluster_num << endl;
            Q.push(*iter);
            (*iter)->queued = true;
            while (!(Q.empty()))
            {
                // dequeue
                Center *head_center = Q.front();
                head_center->clusterID = cluster_num;
                for (auto iter : head_center->members)
                {
                    iter->clusterID = cluster_num;
                }
                Q.pop();
                // find reachable neighbor in S*, add them to Queue
                for (auto Ap_iter = head_center->center->Ap.begin(); Ap_iter != head_center->center->Ap.end(); Ap_iter++)
                {

                    if ((*Ap_iter)->members.size() >= MinPts)
                    {
                        if ((*Ap_iter)->queued == false)
                        {

                            if (BCP_edit(head_center->members, (*Ap_iter)->members, epsilon, distance_matrix) <= epsilon)
                            {
                                Q.push(*Ap_iter);
                                (*Ap_iter)->queued = true;
                            }
                        }
                    }
                    else
                    {

                        for (auto ApMember_iter = (*Ap_iter)->members.begin(); ApMember_iter != (*Ap_iter)->members.end(); ApMember_iter++)
                        {
                            if ((*ApMember_iter)->is_core == true)
                            {

                                Center *core_center = (*ApMember_iter)->center;
                                if (core_center->queued == false)
                                {
                                    if (BCP_edit(head_center->members, core_center->members, epsilon, distance_matrix) <= epsilon)
                                    {
                                        Q.push(core_center);
                                        core_center->queued = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            cluster_num++;
        }
    }
    cout << "Max cluster number is : " << cluster_num << endl;
    return cluster_num;
}
void cluster_border(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            for (auto C_iter = (*iter)->center->Ap.begin(); C_iter != (*iter)->center->Ap.end(); C_iter++)
            {
                for (auto member_iter = (*C_iter)->members.begin(); member_iter != (*C_iter)->members.end(); member_iter++)
                {
                    if ((*member_iter)->is_core == true)
                    {
                        double dist = distance_matrix[(*iter)->ID][(*member_iter)->ID];
                        if (dist == -1)
                            dist = NormDist(**iter, **member_iter);
                        if (dist <= epsilon)
                        {
                            (*iter)->clusterID = (*member_iter)->clusterID;
                            break;
                        }
                    }
                }
            }
        }
    }
}
void cluster_border_edit(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            for (auto C_iter = (*iter)->center->Ap.begin(); C_iter != (*iter)->center->Ap.end(); C_iter++)
            {
                for (auto member_iter = (*C_iter)->members.begin(); member_iter != (*C_iter)->members.end(); member_iter++)
                {
                    if ((*member_iter)->is_core == true)
                    {
                        double dist = distance_matrix[(*iter)->ID][(*member_iter)->ID];
                        if (dist == -1)
                            dist = LevenshteinDist(**iter, **member_iter);
                        if (dist <= epsilon)
                        {
                            (*iter)->clusterID = (*member_iter)->clusterID;
                            break;
                        }
                    }
                }
            }
        }
    }
}
void cluster_border_approx(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            for (auto C_iter = C.begin(); C_iter != C.end(); C_iter++)
            {
                double dist = distance_matrix[(*iter)->ID][(*C_iter)->ID];
                if (dist == -1)
                    dist = NormDist(**iter, **C_iter);
                if (dist <= (*C_iter)->r + epsilon)
                {
                    Point *ptest = (*iter);
                    Center *ctest = (*C_iter);
                    (*iter)->clusterID = (*C_iter)->clusterID;
                    break;
                }
            }
        }
    }
}
void cluster_border_approx_edit(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix, vector<Center *> &C)
{
    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            for (auto C_iter = C.begin(); C_iter != C.end(); C_iter++)
            {
                double dist = distance_matrix[(*iter)->ID][(*C_iter)->ID];
                if (dist == -1)
                    dist = LevenshteinDist(**iter, **C_iter);
                if (dist <= (*C_iter)->r + epsilon)
                {
                    Point *ptest = (*iter);
                    Center *ctest = (*C_iter);
                    (*iter)->clusterID = (*C_iter)->clusterID;
                    break;
                }
            }
        }
    }
}
void DBSCAN_exact(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename)
{
    double r = epsilon * rho; // kcenter termination condition
    vector<Point *> data;     // dataset

    // load data
    load_data(data, dim, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    vector<Center *> centers;
    auto Alg_start = chrono::high_resolution_clock::now();

    // Kcenter process
    cout << "*******************K-center procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    Kcenter(data, centers, r, distance_matrix); // kcenter过程

    cout << "Kcenter complete! "
         << "Centers number: " << centers.size() << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Kcenter time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Calculate neighbor relationship between centers
    cout << "*******************Calculating neighbors*******************" << endl;
    start = chrono::high_resolution_clock::now();

    calculate_Ap(centers, epsilon, distance_matrix);

    cout << "Calculate Ap complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculate neighbors time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // mark core process
    cout << "*******************Mark core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();
    vector<Center *> C;

    mark_core(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Mark core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Mark core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster core process
    cout << "*******************Cluster core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    int cluster_num = cluster_core_approx(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Cluster core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster border process
    cout << "*******************Cluster border procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    cluster_border(data, epsilon, distance_matrix, C);

    cout << "Cluster border complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Algorithm complete
    auto Alg_end = chrono::high_resolution_clock::now();
    cout << "Complete!" << endl;
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    // statistics
    int outlier_num = 0;
    int core_num = 0;
    int border_num = 0;

    vector<int> clusters(cluster_num, 0);

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            outlier_num++;
        }
        else if ((*iter)->is_core == true)
        {
            core_num++;
            clusters[(*iter)->clusterID]++;
        }
        else
        {
            border_num++;
            clusters[(*iter)->clusterID]++;
        }
    }
    cout << "Core point count: " << core_num << endl;
    cout << "Border point count: " << border_num << endl;
    cout << "Outlier count: " << outlier_num << endl;

    std::vector<int> indices(clusters.size());

    for (int i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              { return clusters[a] > clusters[b]; });

    std::map<int, int> clusterID_map;
    for (int i = 0; i < indices.size(); ++i)
    {
        clusterID_map[indices[i]] = i;
    }

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        auto mapIt = clusterID_map.find((*iter)->clusterID);
        if (mapIt != clusterID_map.end())
        {
            (*iter)->clusterID = mapIt->second;
        }
    }

    sort(clusters.rbegin(), clusters.rend());

    output(data, centers, output_data_filename, output_center_filename);
}
void DBSCAN_exact_covertree(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename, int max_layer, int terminate_layer)
{
    double r = epsilon * rho;                 // kcenter termination condition
    int epsilon_net_layer = floor(log2(r));   // epsilon net的层数
    double r_max = pow(2, epsilon_net_layer); // epsilonnet中每个center的最大半径

    vector<Point *> data; // dataset

    // load data
    load_data(data, dim, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    auto Alg_start = chrono::high_resolution_clock::now();

    // CoverTree建树过程

    cout << "*******************CoverTree construction procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    CoverTree *T = new CoverTree(max_layer, data, distance_matrix);

    cout << "Construction complete! " << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "CoverTree construction time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // 计算epsilon net
    cout << "*******************Calculating epsilon net*******************" << endl;
    start = chrono::high_resolution_clock::now();

    vector<Center *> centers;
    calculate_epsilon_net(T, epsilon_net_layer, distance_matrix, centers);

    cout << "Calculate epsilon net complete! "
         << "Centers number: " << centers.size() << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculating epsilon net time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Calculate neighbor relationship between centers
    cout << "*******************Calculating neighbors*******************" << endl;
    start = chrono::high_resolution_clock::now();

    calculate_Ap_covertree(centers, data, r_max, epsilon, distance_matrix);

    cout << "Calculate Ap complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculate neighbors time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // 至此，centers中存储了epsilon net中各centers的Node，以及每个center的Ap集合。尚未计算每个center的weight，以及具体包含哪些points
    cout << "*******************Mark core*******************" << endl;
    start = chrono::high_resolution_clock::now();
    vector<Center *> C;
    mark_core_covertree(centers, T, distance_matrix, MinPts, epsilon, data, C);
    cout << "Mark core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Mark core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    cout << "*******************connect center*******************" << endl;
    start = chrono::high_resolution_clock::now();
    connect_centers_bruteforce(centers, distance_matrix, data, terminate_layer, r_max, MinPts, epsilon);
    cout << "connect center complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "connect center time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // 简单处理重复点
    for (int i = 0; i < n; i++)
    {
        if (data[i]->clusterID == -1)
        {
            for (int j = 0; j < n; j++)
            {
                if (data[j]->coords == data[i]->coords && i != j)
                {
                    data[i]->clusterID = data[j]->clusterID;
                }
            }
        }
    }

    // 统计过程

    auto Alg_end = chrono::high_resolution_clock::now();
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    cout << centers.size() << " centers in total." << endl;
    int clustered = 0, outlier = 0;

    for (auto item : data)
    {
        if (item->clusterID == -1)
            outlier++;
        else
            clustered++;
    }
    cout << clustered << " points been clustered." << endl;
    cout << outlier << " points are outliers." << endl;
    output(data, centers, output_data_filename, output_center_filename);
}
void DBSCAN_exact_edit(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename)
{
    double r = epsilon * rho; // kcenter termination condition
    vector<Point *> data;     // dataset

    // load data
    load_text_data(data, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    vector<Center *> centers;
    auto Alg_start = chrono::high_resolution_clock::now();

    // Kcenter process
    cout << "*******************K-center procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    Kcenter_edit(data, centers, r, distance_matrix); // kcenter过程

    cout << "Kcenter complete! "
         << "Centers number: " << centers.size() << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Kcenter time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Calculate neighbor relationship between centers
    cout << "*******************Calculating neighbors*******************" << endl;
    start = chrono::high_resolution_clock::now();

    calculate_Ap(centers, epsilon, distance_matrix);

    cout << "Calculate Ap complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculate neighbors time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // mark core process
    cout << "*******************Mark core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();
    vector<Center *> C;

    mark_core_edit(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Mark core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Mark core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster core process
    cout << "*******************Cluster core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    int cluster_num = cluster_core_edit(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Cluster core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster border process
    cout << "*******************Cluster border procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    cluster_border_edit(data, epsilon, distance_matrix, C);

    cout << "Cluster border complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Algorithm complete
    auto Alg_end = chrono::high_resolution_clock::now();
    cout << "Complete!" << endl;
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    // statistics
    int outlier_num = 0;
    int core_num = 0;
    int border_num = 0;

    vector<int> clusters(cluster_num, 0);

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            outlier_num++;
        }
        else if ((*iter)->is_core == true)
        {
            core_num++;
            clusters[(*iter)->clusterID]++;
        }
        else
        {
            border_num++;
            clusters[(*iter)->clusterID]++;
        }
    }
    cout << "Core point count: " << core_num << endl;
    cout << "Border point count: " << border_num << endl;
    cout << "Outlier count: " << outlier_num << endl;

    /* std::vector<int> indices(clusters.size());

    for (int i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              { return clusters[a] > clusters[b]; });

    std::map<int, int> clusterID_map;
    for (int i = 0; i < indices.size(); ++i)
    {
        clusterID_map[indices[i]] = i;
    }

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        auto mapIt = clusterID_map.find((*iter)->clusterID);
        if (mapIt != clusterID_map.end())
        {
            (*iter)->clusterID = mapIt->second;
        }
    }

    sort(clusters.rbegin(), clusters.rend()); */

    output(data, centers, output_data_filename, output_center_filename);
}
void DBSCAN_approx(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename)
{
    double r = epsilon * rho; // kcenter termination condition
    vector<Point *> data;     // dataset

    // load data
    load_data(data, dim, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    vector<Center *> centers;
    auto Alg_start = chrono::high_resolution_clock::now();

    // Kcenter process
    cout << "*******************K-center procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    Kcenter(data, centers, r, distance_matrix); // kcenter过程

    cout << "Kcenter complete! "
         << "Centers number: " << centers.size() << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Kcenter time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Calculate neighbor relationship between centers
    cout << "*******************Calculating neighbors*******************" << endl;
    start = chrono::high_resolution_clock::now();

    calculate_Ap(centers, epsilon, distance_matrix);

    cout << "Calculate Ap complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculate neighbors time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // mark core process
    cout << "*******************Mark core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();
    vector<Center *> C;

    mark_core(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Mark core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Mark core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster core process
    cout << "*******************Cluster core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    int cluster_num = cluster_core_approx(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Cluster core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster border process
    cout << "*******************Cluster border procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    cluster_border_approx(data, epsilon, distance_matrix, C);

    cout << "Cluster border complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Algorithm complete
    auto Alg_end = chrono::high_resolution_clock::now();
    cout << "Complete!" << endl;
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    // statistics
    int outlier_num = 0;
    int core_num = 0;
    int border_num = 0;

    vector<int> clusters(cluster_num, 0);

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            outlier_num++;
        }
        else if ((*iter)->is_core == true)
        {
            core_num++;
            clusters[(*iter)->clusterID]++;
        }
        else
        {
            border_num++;
            clusters[(*iter)->clusterID]++;
        }
    }
    cout << "Core point count: " << core_num << endl;
    cout << "Border point count: " << border_num << endl;
    cout << "Outlier count: " << outlier_num << endl;

    std::vector<int> indices(clusters.size());

    for (int i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              { return clusters[a] > clusters[b]; });

    std::map<int, int> clusterID_map;
    for (int i = 0; i < indices.size(); ++i)
    {
        clusterID_map[indices[i]] = i;
    }

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        auto mapIt = clusterID_map.find((*iter)->clusterID);
        if (mapIt != clusterID_map.end())
        {
            (*iter)->clusterID = mapIt->second;
        }
    }

    sort(clusters.rbegin(), clusters.rend());

    output(data, centers, output_data_filename, output_center_filename);
}
void DBSCAN_approx_random(int dim, int n, int MinPts, double epsilon, double rho, int init_num, int farthest_num, int select_num, string input_fileName, string output_data_filename, string output_center_filename)
{
    double r = epsilon * rho; // kcenter termination condition
    vector<Point *> data;     // dataset

    // load data
    load_data(data, dim, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    vector<Center *> centers;
    auto Alg_start = chrono::high_resolution_clock::now();

    // Kcenter process
    cout << "*******************K-center procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    Randomized_Kcenter(data, centers, r, distance_matrix, init_num, farthest_num, select_num); // kcenter过程

    cout << "Kcenter complete! "
         << "Centers number: " << centers.size() << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Kcenter time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Calculate neighbor relationship between centers
    cout << "*******************Calculating neighbors*******************" << endl;
    start = chrono::high_resolution_clock::now();

    calculate_Ap(centers, epsilon, distance_matrix);

    cout << "Calculate Ap complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculate neighbors time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // mark core process
    cout << "*******************Mark core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();
    vector<Center *> C;

    mark_core(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Mark core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Mark core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster core process
    cout << "*******************Cluster core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    int cluster_num = cluster_core_approx(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Cluster core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster border process
    cout << "*******************Cluster border procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    cluster_border_approx(data, epsilon, distance_matrix, C);

    cout << "Cluster border complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Algorithm complete
    auto Alg_end = chrono::high_resolution_clock::now();
    cout << "Complete!" << endl;
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    // statistics
    int outlier_num = 0;
    int core_num = 0;
    int border_num = 0;

    vector<int> clusters(cluster_num, 0);

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            outlier_num++;
        }
        else if ((*iter)->is_core == true)
        {
            core_num++;
            clusters[(*iter)->clusterID]++;
        }
        else
        {
            border_num++;
            clusters[(*iter)->clusterID]++;
        }
    }
    cout << "Core point count: " << core_num << endl;
    cout << "Border point count: " << border_num << endl;
    cout << "Outlier count: " << outlier_num << endl;

    std::vector<int> indices(clusters.size());

    for (int i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              { return clusters[a] > clusters[b]; });

    std::map<int, int> clusterID_map;
    for (int i = 0; i < indices.size(); ++i)
    {
        clusterID_map[indices[i]] = i;
    }

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        auto mapIt = clusterID_map.find((*iter)->clusterID);
        if (mapIt != clusterID_map.end())
        {
            (*iter)->clusterID = mapIt->second;
        }
    }

    sort(clusters.rbegin(), clusters.rend());

    output(data, centers, output_data_filename, output_center_filename);
}
void DBSCAN_approx_edit(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename)
{
    double r = epsilon * rho; // kcenter termination condition
    vector<Point *> data;     // dataset

    // load data
    load_text_data(data, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    vector<Center *> centers;
    auto Alg_start = chrono::high_resolution_clock::now();

    // Kcenter process
    cout << "*******************K-center procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    Kcenter_edit(data, centers, r, distance_matrix); // kcenter过程

    cout << "Kcenter complete! "
         << "Centers number: " << centers.size() << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Kcenter time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Calculate neighbor relationship between centers
    cout << "*******************Calculating neighbors*******************" << endl;
    start = chrono::high_resolution_clock::now();

    calculate_Ap(centers, epsilon, distance_matrix);

    cout << "Calculate Ap complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculate neighbors time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // mark core process
    cout << "*******************Mark core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();
    vector<Center *> C;

    mark_core_edit(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Mark core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Mark core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster core process
    cout << "*******************Cluster core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    int cluster_num = cluster_core_approx(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Cluster core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster border process
    cout << "*******************Cluster border procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    cluster_border_approx_edit(data, epsilon, distance_matrix, C);

    cout << "Cluster border complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Algorithm complete
    auto Alg_end = chrono::high_resolution_clock::now();
    cout << "Complete!" << endl;
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    // statistics
    int outlier_num = 0;
    int core_num = 0;
    int border_num = 0;

    vector<int> clusters(cluster_num, 0);

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            outlier_num++;
        }
        else if ((*iter)->is_core == true)
        {
            core_num++;
            clusters[(*iter)->clusterID]++;
        }
        else
        {
            border_num++;
            clusters[(*iter)->clusterID]++;
        }
    }
    cout << "Core point count: " << core_num << endl;
    cout << "Border point count: " << border_num << endl;
    cout << "Outlier count: " << outlier_num << endl;

    std::vector<int> indices(clusters.size());

    for (int i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              { return clusters[a] > clusters[b]; });

    std::map<int, int> clusterID_map;
    for (int i = 0; i < indices.size(); ++i)
    {
        clusterID_map[indices[i]] = i;
    }

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        auto mapIt = clusterID_map.find((*iter)->clusterID);
        if (mapIt != clusterID_map.end())
        {
            (*iter)->clusterID = mapIt->second;
        }
    }

    sort(clusters.rbegin(), clusters.rend());

    output(data, centers, output_data_filename, output_center_filename);
}
void DBSCAN_vanilla(int dim, int n, int MinPts, double epsilon, string input_fileName, string output_data_filename, string output_center_filename)
{
    vector<Point *> data; // dataset

    // load data
    load_data(data, dim, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    auto Alg_start = chrono::high_resolution_clock::now();
    // range search process
    cout << "*******************Range search procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    range_search(data, MinPts, epsilon, distance_matrix); // 每个点遍历所有点进行range search,并确定每个点是否为core point

    cout << "Range search complete!" << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Range search time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster core process
    cout << "*******************Clustering procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    int cluster_num = cluster_core_vanilla(data, epsilon, distance_matrix);

    cout << "Clustering complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Clustering time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Algorithm complete
    auto Alg_end = chrono::high_resolution_clock::now();
    cout << "Complete!" << endl;
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    // statistics
    int outlier_num = 0;
    int core_num = 0;
    int border_num = 0;

    vector<int> clusters(cluster_num, 0);

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            outlier_num++;
        }
        else if ((*iter)->is_core == true)
        {
            core_num++;
            clusters[(*iter)->clusterID]++;
        }
        else
        {
            border_num++;
            clusters[(*iter)->clusterID]++;
        }
    }
    cout << "Core point count: " << core_num << endl;
    cout << "Border point count: " << border_num << endl;
    cout << "Outlier count: " << outlier_num << endl;

    std::vector<int> indices(clusters.size());

    for (int i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              { return clusters[a] > clusters[b]; });

    std::map<int, int> clusterID_map;
    for (int i = 0; i < indices.size(); ++i)
    {
        clusterID_map[indices[i]] = i;
    }

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        auto mapIt = clusterID_map.find((*iter)->clusterID);
        if (mapIt != clusterID_map.end())
        {
            (*iter)->clusterID = mapIt->second;
        }
    }

    sort(clusters.rbegin(), clusters.rend());

    vector<Center *> centers;
    output(data, centers, output_data_filename, output_center_filename);
}
void DBSCAN_vanilla_edit(int dim, int n, int MinPts, double epsilon, string input_fileName, string output_data_filename, string output_center_filename)
{
    vector<Point *> data; // dataset

    // load data
    load_text_data(data, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    auto Alg_start = chrono::high_resolution_clock::now();
    // range search process
    cout << "*******************Range search procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    range_search_edit(data, MinPts, epsilon, distance_matrix); // 每个点遍历所有点进行range search,并确定每个点是否为core point

    cout << "Range search complete!" << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Range search time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster core process
    cout << "*******************Clustering procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    int cluster_num = cluster_core_vanilla(data, epsilon, distance_matrix);

    cout << "Clustering complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Clustering time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Algorithm complete
    auto Alg_end = chrono::high_resolution_clock::now();
    cout << "Complete!" << endl;
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    // statistics
    int outlier_num = 0;
    int core_num = 0;
    int border_num = 0;

    vector<int> clusters(cluster_num, 0);

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            outlier_num++;
        }
        else if ((*iter)->is_core == true)
        {
            core_num++;
            clusters[(*iter)->clusterID]++;
        }
        else
        {
            border_num++;
            clusters[(*iter)->clusterID]++;
        }
    }
    cout << "Core point count: " << core_num << endl;
    cout << "Border point count: " << border_num << endl;
    cout << "Outlier count: " << outlier_num << endl;

    std::vector<int> indices(clusters.size());

    for (int i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              { return clusters[a] > clusters[b]; });

    std::map<int, int> clusterID_map;
    for (int i = 0; i < indices.size(); ++i)
    {
        clusterID_map[indices[i]] = i;
    }

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        auto mapIt = clusterID_map.find((*iter)->clusterID);
        if (mapIt != clusterID_map.end())
        {
            (*iter)->clusterID = mapIt->second;
        }
    }

    sort(clusters.rbegin(), clusters.rend());

    vector<Center *> centers;
    output(data, centers, output_data_filename, output_center_filename);
}

void DBSCAN_approx_covertree(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename, int max_layer, int terminate_layer)
{
    double r = epsilon * rho;                 // kcenter termination condition
    int epsilon_net_layer = floor(log2(r));   // epsilon net的层数
    double r_max = pow(2, epsilon_net_layer); // epsilonnet中每个center的最大半径

    vector<Point *> data; // dataset

    // load data
    load_data(data, dim, n, input_fileName);
    cout << "load data complete!" << endl;

    vector<vector<double>> distance_matrix(n, vector<double>(n, -1));
    for (int i = 0; i < n; i++)
    {
        distance_matrix[i][i] = 0;
    }

    auto Alg_start = chrono::high_resolution_clock::now();

    // CoverTree建树过程

    cout << "*******************CoverTree construction procedure*******************" << endl;
    auto start = chrono::high_resolution_clock::now();

    CoverTree *T = new CoverTree(max_layer, data, distance_matrix);

    cout << "Construction complete! " << endl;
    auto end = chrono::high_resolution_clock::now();
    std::cout << "CoverTree construction time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // 计算epsilon net
    cout << "*******************Calculating epsilon net*******************" << endl;
    start = chrono::high_resolution_clock::now();

    vector<Center *> centers;
    calculate_epsilon_net(T, epsilon_net_layer, distance_matrix, centers);

    cout << "Calculate epsilon net complete! "
         << "Centers number: " << centers.size() << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculating epsilon net time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Calculate neighbor relationship between centers
    cout << "*******************Calculating neighbors*******************" << endl;
    start = chrono::high_resolution_clock::now();

    calculate_Ap_covertree(centers, data, r_max, epsilon, distance_matrix);

    cout << "Calculate Ap complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Calculate neighbors time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // 至此，centers中存储了epsilon net中各centers的Node，以及每个center的Ap集合。尚未计算每个center的weight，以及具体包含哪些points
    cout << "*******************Mark core*******************" << endl;
    start = chrono::high_resolution_clock::now();
    vector<Center *> C;
    mark_core_covertree(centers, T, distance_matrix, MinPts, epsilon, data, C);
    cout << "Mark core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Mark core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster core process
    cout << "*******************Cluster core procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    int cluster_num = cluster_core_approx(data, centers, epsilon, MinPts, distance_matrix, C);

    cout << "Cluster core complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster core time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // cluster border process
    cout << "*******************Cluster border procedure*******************" << endl;
    start = chrono::high_resolution_clock::now();

    cluster_border_approx(data, epsilon, distance_matrix, C);

    cout << "Cluster border complete!" << endl;
    end = chrono::high_resolution_clock::now();
    std::cout << "Cluster border time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
    cout << "********************************************************" << endl;

    // Algorithm complete
    auto Alg_end = chrono::high_resolution_clock::now();
    cout << "Complete!" << endl;
    std::cout << "Total time: " << std::chrono::duration<double>(Alg_end - Alg_start).count() * 1000 << " ms" << std::endl;

    // statistics
    int outlier_num = 0;
    int core_num = 0;
    int border_num = 0;

    vector<int> clusters(cluster_num, 0);

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        if ((*iter)->clusterID == -1)
        {
            outlier_num++;
        }
        else if ((*iter)->is_core == true)
        {
            core_num++;
            clusters[(*iter)->clusterID]++;
        }
        else
        {
            border_num++;
            clusters[(*iter)->clusterID]++;
        }
    }
    cout << "Core point count: " << core_num << endl;
    cout << "Border point count: " << border_num << endl;
    cout << "Outlier count: " << outlier_num << endl;

    std::vector<int> indices(clusters.size());

    for (int i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              { return clusters[a] > clusters[b]; });

    std::map<int, int> clusterID_map;
    for (int i = 0; i < indices.size(); ++i)
    {
        clusterID_map[indices[i]] = i;
    }

    for (auto iter = data.begin(); iter != data.end(); iter++)
    {
        auto mapIt = clusterID_map.find((*iter)->clusterID);
        if (mapIt != clusterID_map.end())
        {
            (*iter)->clusterID = mapIt->second;
        }
    }

    sort(clusters.rbegin(), clusters.rend());

    output(data, centers, output_data_filename, output_center_filename);
}