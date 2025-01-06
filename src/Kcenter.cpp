
#include "../header/Kcenter.h"
using namespace std;

bool cmp(dist_struct a, dist_struct b)
{
    return *(a.dists_ptr) > *(b.dists_ptr);
}
void Randomized_Kcenter(vector<Point *> &data, vector<Center *> &centers, double r, vector<vector<double>> &distance_matrix, int init_num, int farthest_num, int select_num)
{


    int n = data.size();
    srand(time(0));
    int init_center_ID;

    vector<Center *> temp_center(n);   
    vector<double> dists(n, MAX_DIST); 

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


    vector<Point *> candidate_centers; 
    vector<Center *> new_centers;      
    vector<int> shuffle_n;             
    vector<int> shuffle_farthest_num;  
    for (int i = 0; i < n; i++)
        shuffle_n.push_back(i);
    for (int i = 0; i < farthest_num; i++)
        shuffle_farthest_num.push_back(i);
    random_shuffle(shuffle_n.begin(), shuffle_n.end());
    for (int i = 0; i < init_num; i++)
    {
        init_center_ID = shuffle_n[i]; 
        dists[init_center_ID] = 0;
        candidate_centers.push_back(data[init_center_ID]);
    }



    double farthest_set_to_E = MAX_DIST; 
    do
    {
        
        for (auto candidate : candidate_centers)
        {
            Center *new_center = new Center(*candidate);
            candidate->center = new_center;
            new_center->center = new_center;
            new_center->members.push_back(candidate); 
            new_centers.push_back(new_center);
            centers.push_back(new_center);
        }



        // farthest_set_to_E = MAX_DIST;
        //  auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < n; i++)
        {
            if (data[i]->center == NULL)
            { 

                // auto start = chrono::high_resolution_clock::now();
                double dist_to_new_centers = MAX_DIST; 
                for (auto center_it : new_centers)
                {
                    double dist_to_center = NormDist(*data[i], *center_it);
                    distance_matrix[i][center_it->ID] = dist_to_center;
                    distance_matrix[center_it->ID][i] = dist_to_center; 
                    if (dist_to_center < dist_to_new_centers)

                    { 
                        dist_to_new_centers = dist_to_center;
                        if (dist_to_center < dists[i])
                            temp_center[i] = center_it;
                    }
                }

                if (dist_to_new_centers < dists[i])
                    dists[i] = dist_to_new_centers;
            }
        }

        
        candidate_centers.clear();
        new_centers.clear();

        // sort(dists_ptrs.begin(), dists_ptrs.end(), cmp);
        nth_element(dists_ptrs.begin(), dists_ptrs.begin() + farthest_num, dists_ptrs.end(), cmp);
        
        if (*(dists_ptrs[farthest_num - 1].dists_ptr) < farthest_set_to_E)
            farthest_set_to_E = *(dists_ptrs[farthest_num - 1].dists_ptr);

        
        random_shuffle(shuffle_farthest_num.begin(), shuffle_farthest_num.end());

        for (int i = 0; i < select_num; i++)
        {
            init_center_ID = dists_ptrs[shuffle_farthest_num[i]].ID;
            dists[init_center_ID] = 0;
            candidate_centers.push_back(data[init_center_ID]);
        }

        loop++;

        

        // end = chrono::high_resolution_clock::now();
        // std::cout << "Origin method total loop time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

    } while (farthest_set_to_E > r && loop <= n);
    
    for (int i = 0; i < farthest_num; i++)
    {
        Center *new_center = new Center(*(data[dists_ptrs[i].ID]));
        data[dists_ptrs[i].ID]->center = new_center;
        new_center->center = new_center;
        new_center->members.push_back(data[dists_ptrs[i].ID]); 
        centers.push_back(new_center);
    }

    
    for (int i = 0; i < n; i++)
    {
        if (data[i]->center == NULL)
        { 

            data[i]->center = temp_center[i];

            
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
    int count = 0;
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
                // distance_matrix[i][new_center->ID] = dist_to_new_center;
                // distance_matrix[new_center->ID][i] = dist_to_new_center;
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
        count++;
        if (count >= 100)
        {
            count = 0;
            cout << loop << " loops solved." << endl;
            cout << "max_r = " << max_r << "; r = " << r << endl;
        }
        // cout << "center count: " << loop << endl;
        // cout << "max_r = " << max_r << "; r = " << r << endl;

        // end = chrono::high_resolution_clock::now();
        // std::cout << "Origin method total loop time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

    } while (max_r > r && loop <= n);

    // assign remaining points to their closest center
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

void Kcenter_angular(vector<Point *> &data, vector<Center *> &centers, double r, vector<vector<double>> &distance_matrix)
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
    int count = 0;
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
                double dist_to_new_center = CosinSim(*data[i], *new_center);
                // auto end = chrono::high_resolution_clock::now();
                // std::cout << "count one dist time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

                // start = chrono::high_resolution_clock::now();
                // distance_matrix[i][new_center->ID] = dist_to_new_center;
                // distance_matrix[new_center->ID][i] = dist_to_new_center;
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
        count++;
        if (count >= 100)
        {
            count = 0;
            cout << loop << " loops solved." << endl;
            cout << "max_r = " << max_r << "; r = " << r << endl;
        }
        // cout << "center count: " << loop << endl;
        // cout << "max_r = " << max_r << "; r = " << r << endl;

        // end = chrono::high_resolution_clock::now();
        // std::cout << "Origin method total loop time: " << std::chrono::duration<double>(end - start).count() * 1000 << " ms" << std::endl;

    } while (max_r > r && loop <= n);

    // assign remaining points to their closest center
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
