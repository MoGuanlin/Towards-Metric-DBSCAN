
#include "../header/Kcenter.h"
using namespace std;

bool cmp(dist_struct a, dist_struct b)
{
    return *(a.dists_ptr) > *(b.dists_ptr);
}
void Randomized_Kcenter(vector<Point *> &data, vector<Center *> &centers, double r, vector<vector<double>> &distance_matrix, int init_num, int farthest_num, int select_num)
{

    // 随机选择init_num个点作为初始的center
    int n = data.size();
    srand(time(0));
    int init_center_ID;

    vector<Center *> temp_center(n);   // 存储每个点所属临时center
    vector<double> dists(n, MAX_DIST); // 存储每个点与centers集合的距离,初始化为MAX_DIST

    // dist_ptrs存储double *,ID结构体
    vector<dist_struct> dists_ptrs;
    for (int i = 0; i < n; i++)
    {
        dist_struct t;
        t.ID = i;
        t.dists_ptr = &dists[i];
        dists_ptrs.push_back(t);
    }

    int loop = 1;
    double max_r = 0;

    // 把初始center加入center集合(将data中的相应数据点的center属性赋正值，作为与普通point的区别)
    vector<Point *> candidate_centers; // point形式的新center
    vector<Center *> new_centers;      // center形式的新center
    vector<int> shuffle_n;             // 用于在n个数里面取init_num个互不相同的随机数
    vector<int> shuffle_farthest_num;  // 用于在farthest_num个数里面取select_num个互不相同的随机数
    for (int i = 0; i < n; i++)
        shuffle_n.push_back(i);
    for (int i = 0; i < farthest_num; i++)
        shuffle_farthest_num.push_back(i);
    random_shuffle(shuffle_n.begin(), shuffle_n.end());
    for (int i = 0; i < init_num; i++)
    {
        init_center_ID = shuffle_n[i]; // 注意避免重复序号，n个里面选init_num个
        dists[init_center_ID] = 0;
        candidate_centers.push_back(data[init_center_ID]);
    }

    // 每一轮迭代中，遍历所有数据点，找出与centers集合中的点最远的那个点
    // 每个数据点临时存储它与centers集合的距离(最短距离)，和这个最短距离对应的center(即这个点的临时所属center)
    // 对于新的center，每个点更新最短距离值和最短距离对应的center

    double farthest_set_to_E = MAX_DIST; // 记录下一批最远的(1+δ)z个点与当前已选center集E的距离
    do
    {
        // 把选出的所有candidate center添加到centers集合中
        for (auto candidate : candidate_centers)
        {
            Center *new_center = new Center(*candidate);
            candidate->center = new_center;
            new_center->center = new_center;
            new_center->members.push_back(candidate); // 向member中添加的应是原始数据点类
            new_centers.push_back(new_center);
            centers.push_back(new_center);
        }

        // 每次要算新center到其它非center点的距离：
        // 大matrix data：784*70000
        // 每次为new_center创建一个向量：784*1，把这个向量的转置与大矩阵进行点乘，
        // 获得一个向量1*70000，这个向量的每一维值就是相应ID数据点与new_center之间的欧氏距离

        // 计算new_center与每个点的距离,并更新每个点的与center集距离，以及临时所属center

        // farthest_set_to_E = MAX_DIST;
        //  auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < n; i++)
        {
            if (data[i]->center == NULL)
            { // 仅操作当前不是center的point

                // auto start = chrono::high_resolution_clock::now();
                double dist_to_new_centers = MAX_DIST; // 该点与新的centers之间的最短距离
                for (auto center_it : new_centers)
                {
                    double dist_to_center = NormDist(*data[i], *center_it);
                    distance_matrix[i][center_it->ID] = dist_to_center;
                    distance_matrix[center_it->ID][i] = dist_to_center; // 更新距离矩阵
                    if (dist_to_center < dist_to_new_centers)

                    { // 更新新一批centers与点i的最近距离
                        dist_to_new_centers = dist_to_center;
                        if (dist_to_center < dists[i])
                            // 某个新center既更新了当前批center与点i的距离，又与点i的距离更近于所有的老center，则更新i的所属center
                            temp_center[i] = center_it;
                    }
                }
                // 根据新一批centers与点i的最近距离，更新dist[i]
                if (dist_to_new_centers < dists[i])
                    dists[i] = dist_to_new_centers;
            }
        }

        // 更新完所有点当前与center集E的距离，先清空candidate centers与new centers,准备下一轮重新添加
        candidate_centers.clear();
        new_centers.clear();

        // 接下来要找到(1+δ)z个(farthest num个)距离E最远的点
        // 然后在这farthest num个数中随机选select num个ID，把这些ID对应的点加入candidate centers，进入下一轮循环
        // sort(dists_ptrs.begin(), dists_ptrs.end(), cmp);
        nth_element(dists_ptrs.begin(), dists_ptrs.begin() + farthest_num, dists_ptrs.end(), cmp);
        // 根据sort结果，更新farthest_set_to E
        if (*(dists_ptrs[farthest_num - 1].dists_ptr) < farthest_set_to_E)
            farthest_set_to_E = *(dists_ptrs[farthest_num - 1].dists_ptr);

        // 打乱，选前farthest num中随机的select num个，成为新的center
        random_shuffle(shuffle_farthest_num.begin(), shuffle_farthest_num.end());

        for (int i = 0; i < select_num; i++)
        {
            init_center_ID = dists_ptrs[shuffle_farthest_num[i]].ID;
            dists[init_center_ID] = 0;
            candidate_centers.push_back(data[init_center_ID]);
        }

        loop++;

        // 迭代完成之后，candidate_center就是新一轮的new_center,将它升级,并加入centers集合中

        // end = chrono::high_resolution_clock::now();
        // std::cout << "Origin method total loop time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

    } while (farthest_set_to_E > r && loop <= n);
    // TODO:将当前dists_ptrs中前farest num个点作为outlier加入centers中，并更新每个point的center属性，统计weight和memberID
    //
    // 将当前dists_ptrs中前farest num个点作为outlier加入centers中

    for (int i = 0; i < farthest_num; i++)
    {
        Center *new_center = new Center(*(data[dists_ptrs[i].ID]));
        data[dists_ptrs[i].ID]->center = new_center;
        new_center->center = new_center;
        new_center->members.push_back(data[dists_ptrs[i].ID]); // 向member中添加的应是原始数据点类
        centers.push_back(new_center);
    }

    //  把temp center设为每个point的center
    for (int i = 0; i < n; i++)
    {
        if (data[i]->center == NULL)
        { // 仅操作当前不是center的point

            data[i]->center = temp_center[i];

            // 更新对应center的members,weight,r
            data[i]->center->members.push_back(data[i]);
            data[i]->center->weight++;
            if (dists[i] > data[i]->center->r)
                data[i]->center->r = dists[i];
        }
    }
}
void Kcenter(vector<Point *> &data, vector<Center *> &centers, double r, vector<vector<double>> &distance_matrix)
{

    // init center
    int n = data.size();
    srand(time(0));
    int init_center_num = rand() % n;
    vector<Center *> temp_center(n);
    vector<double> dists(n, MAX_DIST);
    int loop = 1;
    double max_r = 0;

    Point *candidate_center;
    Center *new_center;
    candidate_center = data[init_center_num];

    do
    {

        new_center = new Center(*candidate_center);
        candidate_center->center = new_center;
        new_center->center = new_center;
        new_center->members.push_back(candidate_center);
        centers.push_back(new_center);

        max_r = 0;
        // auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < n; i++)
        {
            if (data[i]->center == NULL)
            {

                // auto start = chrono::high_resolution_clock::now();
                double dist_to_new_center = NormDist(*data[i], *new_center);
                // auto end = chrono::high_resolution_clock::now();
                // std::cout << "count one dist time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

                // start = chrono::high_resolution_clock::now();
                distance_matrix[i][new_center->ID] = dist_to_new_center;
                distance_matrix[new_center->ID][i] = dist_to_new_center;
                // end = chrono::high_resolution_clock::now();
                // std::cout << "save to matrix time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

                if (dist_to_new_center < dists[i])
                { // renew i stats
                    temp_center[i] = new_center;
                    dists[i] = dist_to_new_center;
                }
                if (dists[i] > max_r)
                { // i is new candidate
                    max_r = dists[i];
                    candidate_center = data[i];
                }
            }
        }
        // auto end = chrono::high_resolution_clock::now();
        // std::cout << "Origin method time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
        loop++;

        // end = chrono::high_resolution_clock::now();
        // std::cout << "Origin method total loop time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

    } while (max_r > r && loop <= n);

    // new center to E
    for (int i = 0; i < n; i++)
    {
        if (data[i]->center == NULL)
        {

            data[i]->center = temp_center[i];

            // update members,weight,r
            data[i]->center->members.push_back(data[i]);
            data[i]->center->weight++;
            if (dists[i] > data[i]->center->r)
                data[i]->center->r = dists[i];
        }
    }
}

void Kcenter_edit(vector<Point *> &data, vector<Center *> &centers, double r, vector<vector<double>> &distance_matrix)
{

    // init center
    int n = data.size();
    srand(time(0));
    int init_center_num = rand() % n;
    vector<Center *> temp_center(n);
    vector<double> dists(n, MAX_DIST);
    int loop = 1;
    double max_r = 0;

    Point *candidate_center;
    Center *new_center;
    candidate_center = data[init_center_num];

    do
    {
        cout << "New center found! center number: " << loop << endl;
        new_center = new Center(*candidate_center); //////////////////
        candidate_center->center = new_center;
        new_center->center = new_center;
        new_center->members.push_back(candidate_center);
        centers.push_back(new_center);

        max_r = 0;
        // auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < n; i++)
        {
            if (data[i]->center == NULL)
            {

                // auto start = chrono::high_resolution_clock::now();
                double dist_to_new_center = LevenshteinDist(*data[i], *new_center);
                // auto end = chrono::high_resolution_clock::now();
                // std::cout << "count one dist time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

                // start = chrono::high_resolution_clock::now();
                distance_matrix[i][new_center->ID] = dist_to_new_center;
                distance_matrix[new_center->ID][i] = dist_to_new_center;
                // end = chrono::high_resolution_clock::now();
                // std::cout << "save to matrix time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

                if (dist_to_new_center < dists[i])
                { // renew i stats
                    temp_center[i] = new_center;
                    dists[i] = dist_to_new_center;
                }
                if (dists[i] > max_r)
                { // i is new candidate
                    max_r = dists[i];
                    candidate_center = data[i];
                }
            }
        }
        // auto end = chrono::high_resolution_clock::now();
        // std::cout << "Origin method time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;
        loop++;
        cout << "max r: " << max_r << endl;

        // end = chrono::high_resolution_clock::now();
        // std::cout << "Origin method total loop time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

    } while (max_r > r && loop <= n);

    // new center to E
    for (int i = 0; i < n; i++)
    {
        if (data[i]->center == NULL)
        {

            data[i]->center = temp_center[i];

            // update members,weight,r
            data[i]->center->members.push_back(data[i]);
            data[i]->center->weight++;
            if (dists[i] > data[i]->center->r)
                data[i]->center->r = dists[i];
        }
    }
}
