// g++ -O3 -w -std=c++11 ./src/*.cpp -o DBSCAN
// g++ -O3 -w -std=c++11 ./src/*.cpp -o DBSCAN -DWIN
//./dbscan "testfile2.txt" "result.txt" 11837 498 10 80 75 100 15
//./dbscan -algo 1 -n 11837 -d 2 -r 10 -k 80 -ds "testfile.txt" -rf "result.txt"

#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <memory.h>
#include <stdio.h>
#include <fstream>
#include <queue>
#include <algorithm>

using namespace std;

// clock_t tt;
// clock_t start;
int DIM;
int cata;

class Point
{
public:
	int *coords; // 坐标
	int clustringID;
	Point **neighbor;
	int countn;

	bool fromE; // 属于P，不属于E，则不是core，返回false
	bool isCorePoint;
	bool visited;

	Point(int *coords)
	{
		this->coords = coords;
		this->fromE = false;
		this->isCorePoint = false;
		this->visited = false;
		clustringID = -1;
		neighbor = NULL;
		countn = 0;
	}
};

bool comp(const Point *p1, const Point *p2)
{

	const int *v1 = p1->coords;
	const int *v2 = p2->coords;
	int pos = 0;
	int curMax = 0;
	int temp = 0;
	int d = DIM;

	for (int i = 0; i < d; i++)
	{
		temp = v1[i] ^ v2[i];
		if (curMax < temp && curMax < (temp ^ curMax))
		{
			// Record the dimension where the different most significant bit lies.
			pos = i;
			curMax = temp;
		}
	}
	return v1[pos] < v2[pos];
}

class Rectangle
{

public:
	// The list of lower coordinates.
	int *minValues;

	// The list of higher coordinates.
	int *maxValues;

	Rectangle()
	{
		this->minValues = new int[DIM];
		this->maxValues = new int[DIM];
	}

	Rectangle(int *minValues, int *maxValues)
	{
		this->minValues = minValues;
		this->maxValues = maxValues;
	}

	~Rectangle()
	{
		delete[] this->maxValues;
		delete[] this->minValues;
	}

	Rectangle &enlarge(const Rectangle &rec)
	{
		int d = DIM;
		for (int i = 0; i < d; i++)
		{
			//(this->minValues[i] < rec.minValues[i]) ?  : this->minValues[i] = rec.minValues[i];
			//(this->maxValues[i] > rec.maxValues[i]) ?  : this->maxValues[i] = rec.maxValues[i];
			// this->minValues[i] = min(this->minValues[i], rec.minValues[i]);
			// this->maxValues[i] = max(this->maxValues[i], rec.maxValues[i]);
			(this->minValues[i] < rec.minValues[i]) ? this->minValues[i] = this->minValues[i] : this->minValues[i] = rec.minValues[i];
			(this->maxValues[i] > rec.maxValues[i]) ? this->maxValues[i] = this->maxValues[i] : this->maxValues[i] = rec.maxValues[i];
		}
		return *(this);
	}

	void setMinCoords(const int *minCoords)
	{
		int d = DIM;
		if (this->minValues == NULL)
		{
			this->minValues = new int[DIM];
		}
		for (int i = 0; i < d; i++)
		{
			this->minValues[i] = minCoords[i];
		}
	}

	void setMaxCoords(const int *maxCoords)
	{
		int d = DIM;
		if (this->maxValues == NULL)
		{
			this->maxValues = new int[DIM];
		}
		for (int i = 0; i < d; i++)
		{
			this->maxValues[i] = maxCoords[i];
		}
	}

	int stateWithSphere(Point *q, double sqr_eps)
	{
		int d = DIM;
		unsigned int closestDist = 0;
		unsigned int farthestDist = 0;
		int temp = 0;
		int temp2 = 0;
		unsigned int sqr_temp = 0;
		unsigned int sqr_temp2 = 0;

		// Find the distances from the closest and farthest points to q in this grid cell.
		for (int i = 0; i < d; i++)
		{
			temp = this->minValues[i] - q->coords[i];
			temp2 = this->maxValues[i] - q->coords[i];
			sqr_temp = temp * temp;
			sqr_temp2 = temp2 * temp2;

			if (temp > 0)
			{
				// q is to the left of this rectangle in this dimension
				closestDist += sqr_temp;
			}
			else if (temp2 < 0)
			{
				// q is to the right of this rectangle in this dimension
				closestDist += sqr_temp2;
			}
			farthestDist += (sqr_temp <= sqr_temp2 ? sqr_temp2 : sqr_temp);
		}

		if (closestDist <= sqr_eps)
		{
			if (farthestDist <= sqr_eps)
				return 1; // fully inside
			return 0;	  // intersect
		}
		return -1; // fully outside
	}
};

class RtreeNode
{
protected:
	//  The list of child nodes.
	vector<RtreeNode *> childNodes;

public:
	// The level of this node in the R-tree. All the leaves are at level-0.
	int level;

	//  The minimum bounding rectangle of the child nodes.
	Rectangle mbr;

	// The data object in of this node if this is a leaf.
	Point *pt;

	RtreeNode()
	{
		// Means that this node currently is empty and invalid.
		this->level = -1;
		this->pt = NULL;
	}

	RtreeNode(Point *pt)
	{
		this->level = 0;
		this->pt = pt;
		this->mbr.setMinCoords(pt->coords);
		this->mbr.setMaxCoords(pt->coords);
	}

	~RtreeNode()
	{
		this->releaseSpace();
	}

	void releaseSpace()
	{
		int num = this->childNodes.size();
		if (num == 0)
			return;
		for (int i = 0; i < num; i++)
		{
			this->childNodes[i]->releaseSpace();
			delete (this->childNodes[i]);
			this->childNodes[i] = NULL;
		}
		vector<RtreeNode *>().swap(this->childNodes);
		this->childNodes.clear();
	}

	vector<RtreeNode *> &getChildNodes()
	{
		return this->childNodes;
	}

	void addChildNode(RtreeNode *child)
	{
		// Update the MBR for the parent node.
		if (this->level == -1)
		{
			this->mbr.setMinCoords(child->mbr.minValues);
			this->mbr.setMaxCoords(child->mbr.maxValues);
			this->level = child->level + 1;
		}

		this->mbr.enlarge(child->mbr);
		this->childNodes.emplace_back(child);
	}
};

class Rtree
{
protected:
	RtreeNode *root;

public:
	Rtree()
	{
		this->root = NULL;
	}

	~Rtree()
	{
		this->releaseSpace();
	}

	void releaseSpace()
	{
		if (this->root != NULL)
		{
			// Release the space of R-tree.
			this->root->releaseSpace();
			delete (this->root);
			this->root = NULL;
		}
	}

	void bulkLoadRtree(vector<Point *> &ptList)
	{
		int fanout = 10;
		int num = ptList.size();
		if (num == 0)
		{
			// showError("Error in Rtree::bulkLoadRtree: the point list is empty!\n");
			return;
		}
		else if (num == 1)
		{
			this->root = new RtreeNode(ptList[0]);
			return;
		}

		queue<RtreeNode *> que, que2;
		queue<RtreeNode *> &work = que;
		queue<RtreeNode *> &store = que2;

		/*************** Sort directly by the coordinates of points ****************/

		sort(ptList.begin(), ptList.end(), comp);
		for (int i = 0; i < num; i++)
		{
			que.push(new RtreeNode(ptList[i]));
		}

		RtreeNode *parent = NULL;
		bool flag = false;
		while (work.size() + store.size() != 1)
		{
			parent = new RtreeNode();

			for (int i = 0; i < fanout; i++)
			{
				if (!work.empty())
				{
					parent->addChildNode(work.front());
					work.pop();
				}
				else
				{
					flag = true;
					break;
				}
			}

			if (flag)
			{
				if (parent->getChildNodes().size() > 0)
				{
					store.push(parent);
				}
				else
				{
					delete (parent);
				}
				std::swap(work, store);
				flag = false;
			}
			else
			{
				store.push(parent);
			}
		}
		if (work.empty())
			this->root = store.front();
		else
			this->root = work.front();
	}

	void rangeQuery(Point *pt, double sqr_eps, vector<Point *> &targetPlace)
	{
		/*if (this->root == NULL) {
			//showError("Error in Rtree::rangeQuery: The root is NULL!\n");
		}*/

		int d = DIM;
		register unsigned int temp;
		// double temp;
		register unsigned int dist = 0;
		Point *pt2;

		// Clear the targetPlace.
		targetPlace.clear();

		int childNum = 0;
		int state = 0;
		unsigned int sqrDist = 0;

		queue<RtreeNode *> que;
		que.push(this->root);
		RtreeNode *cur_node = NULL;

		while (!que.empty())
		{
			cur_node = que.front();
			que.pop();

			// If the current node is a leaf, compute the distance from query point to the point inside this leaf.
			if (cur_node->level == 0)
			{
				// sqrDist = SqrDist(pt ,cur_node->getPoint() );
				pt2 = cur_node->pt;
				dist = 0;
				for (int i = 0; i < d; ++i)
				{
					temp = pt->coords[i] - pt2->coords[i];
					dist += temp * temp;
				}
				sqrDist = dist;

				if (sqrDist <= sqr_eps)
					targetPlace.emplace_back(cur_node->pt);
			}
			else
			{
				// If the current node is not a leaf, check its child nodes and add those intersecting with the rectangle of q to the queue.
				vector<RtreeNode *> &childList = cur_node->getChildNodes();
				childNum = childList.size();
				for (int i = 0; i < childNum; i++)
				{
					//				start1 = getCurrentTime();
					state = childList[i]->mbr.stateWithSphere(pt, sqr_eps);
					//				end1 = getCurrentTime();

					if (state != -1)
					{
						// intersect or fully inside, then add this tree node to queue
						que.push(childList[i]);
					}
				}
			}
		}
	}
};

void readDatasetFromFile(char *filePath, vector<Point *> &pt, int n, int d);
Point *readPointFromFile(FILE *fileHandle, int d);
int CalE(vector<Point *> &points, vector<int> &E, int *PtoE_index, unsigned int *PtoE_value, double r, int init, int k, int ran_sel, int E_length, int n, int d);
void CorePoint(vector<Point *> &points, int *PtoE_index, unsigned int *PtoE_value, vector<int> &E, int E_length, int MinPts, double Eps, double r, int n, int d);
void DBSCAN(vector<Point *> &points, int n);
void AddtoE(vector<Point *> &pt, vector<int> &E, int *PtoE_index, unsigned int *PtoE_value, int *E_new, int E_length, int act_sel);
void init_PtoE(vector<Point *> &pt, int *PtoE_index, unsigned int *PtoE_value, int *E_new, int E_length, int act_sel, int n, int d);
void update_PtoE(vector<Point *> &pt, int *PtoE_index, unsigned int *PtoE_value, int *E_new, int E_length, int act_sel, int n, int d);
unsigned int square_dist(vector<Point *> &pt, int a, int b, int d);
unsigned int furthest_KinP(int *Q, vector<Point *> &pt, unsigned int *PtoE_value, int k, int n);
void AdjustDown(int *top, unsigned int *topvalue, int root, int size);
void select_random_1(int *pick, int num, int sup);
void select_random_2(int *pick, int num, int sup);
int delete_nearE(vector<Point *> &pt, int *E_new, int beforDel, double r, int d);
void generate_distE(unsigned int *dist_E, vector<Point *> &points, vector<int> &E, int E_length, int d);
inline unsigned int query_distE(unsigned int *dist_E, int a, int b, int E_length);
void writeClusterToFile(vector<Point *> &points, char *filePath, int n);

int main(int argc, char *argv[])
{
	// 定义
	int n;
	int d;
	vector<Point *> points;
	vector<int> E;
	int E_length = E.size();
	double r;
	int init;	 // 初次
	int k;		 // 循环第一步
	int ran_sel; // 循环第二步
	double Eps;
	int MinPts;
	///////////////////////////////////////////////////////////////////////////
	/*char data_file[50] = "testfile2.txt";
	char result_file[50] = "result.txt";
	n = 11837;
	d = 498;
	Eps = 10;
	MinPts = 80;
	init = 75;
	k = 100;
	ran_sel = init;
	r = 15;*/

	char *data_file = argv[1];
	char *result_file = argv[2];
	n = atoi(argv[3]);
	d = atoi(argv[4]);
	Eps = atof(argv[5]);
	MinPts = atoi(argv[6]);
	init = atoi(argv[7]);
	k = atoi(argv[8]);
	ran_sel = init;
	r = atof(argv[9]);
	//////////////////////////////////////////////////////////////////////////
	DIM = d;
	int *PtoE_index = new int[n];
	unsigned int *PtoE_value = new unsigned int[n];
	// vector<int> Pz;

	// 执行程序
	// 读入
	readDatasetFromFile(data_file, points, n, d);

	clock_t start = clock();
	//////////////////////////////////////////////////////////////////////////
	// Part-one
	E_length = CalE(points, E, PtoE_index, PtoE_value, r, init, k, ran_sel, E_length, n, d);
	cout << "center num: " << E_length << endl;
	// 释放内存
	// delete[] PtoE_value;
	// PtoE_value = NULL;

	clock_t tt = clock();

	// Part-two-a
	CorePoint(points, PtoE_index, PtoE_value, E, E_length, MinPts, Eps, r, n, d);
	// 释放内存
	PtoE_index = NULL;
	PtoE_value = NULL;
	vector<int>().swap(E);
	// vector<int>().swap(Pz);

	// Part-two-b
	DBSCAN(points, n);
	///////////////////////////////////////////////////////////////////////////
	clock_t end = clock();

	// statistics
	int tocount = 0;
	int *count = new int[cata];
	memset(count, 0, sizeof(int) * cata);
	for (int i = 0; i < points.size(); i++)
	{
		if (points[i]->clustringID != -1)
		{
			count[points[i]->clustringID - 1]++;
			tocount++;
		}
	}
	/*for(int i = 0 ; i < cata ; i++){
		cout << i+1 << "\t" << count[i] << endl;
	}*/
	cout << "Cluster Points\t" << tocount << "/" << n << endl;

	// Print time-consuming
	cout << "part1 time\t" << (double)(tt - start) / (double)CLOCKS_PER_SEC << endl;
	cout << "part2 time\t" << (double)(end - tt) / (double)CLOCKS_PER_SEC << endl;
	cout << "total time\t" << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;

	// 输出
	writeClusterToFile(points, result_file, n);
	// system("pause");

	return 0;
}

// Part 1计算球心
int CalE(vector<Point *> &points, vector<int> &E, int *PtoE_index, unsigned int *PtoE_value, double r, int init, int k, int ran_sel, int E_length, int n, int d)
{

	double sqr_r = r * r;
	double del_r = (r / 2.0) * (r / 2.0);

	// 初次选点
	int *E_first = new int[init]; // 数组值为在Points中的下标
	select_random_1(E_first, init, n);
	// int act_sel = delete_nearE(points, E_first, init, del_r, d);
	int act_sel = init;

	// 加入E，移除P
	AddtoE(points, E, PtoE_index, PtoE_value, E_first, E_length, act_sel);
	// 计算PE距离
	init_PtoE(points, PtoE_index, PtoE_value, E_first, E_length, act_sel, n, d);
	E_length += act_sel;

	// 释放E_first
	delete[] E_first;
	E_first = NULL;

	// 循环
	unsigned int dis_QtoE;
	int *Q = new int[k];								   // 值为在points中的下标
	int *E_new = new int[ran_sel];						   // 值为在points中的下标
	int *index_Q = new int[ran_sel];					   // 值为在Q中的下标
	dis_QtoE = furthest_KinP(Q, points, PtoE_value, k, n); // cout << dis_QtoE << "\t" << E_length << endl;
	while (dis_QtoE > sqr_r)
	{
		select_random_2(index_Q, ran_sel, k);
		for (int i = 0; i < ran_sel; ++i)
		{
			E_new[i] = Q[index_Q[i]];
		}
		act_sel = ran_sel;
		// act_sel = delete_nearE(points, E_new, ran_sel, del_r, d);
		AddtoE(points, E, PtoE_index, PtoE_value, E_new, E_length, act_sel);
		update_PtoE(points, PtoE_index, PtoE_value, E_new, E_length, act_sel, n, d);
		E_length += act_sel;
		dis_QtoE = furthest_KinP(Q, points, PtoE_value, k, n);
		// cout << dis_QtoE <<"\t"<< E_length << endl;
	}
	/*for (int i = 0; i < k; ++i) {
		if(PtoE_value[Q[i]] > sqr_r){
			Pz.emplace_back(Q[i]);
		}
	}*/
	delete[] Q;
	delete[] index_Q;
	delete[] E_new;
	Q = NULL;
	index_Q = NULL;
	E_new = NULL;

	return E_length;
}

// 计算core point

/*void CorePoint(vector<Point*>& points, int* PtoE_index, unsigned int* PtoE_value, vector<int>& E, int E_length, int MinPts, double Eps, double r, int n, int d) {

	bool outliers = false;

	//calculate dist_E
	unsigned int* dist_E = new unsigned int[(E_length - 1) * E_length / 2];
	generate_distE(dist_E, points, E, E_length, d);


	double sqr_r = r * r;
	//calculate EtoP
	vector<int>* EtoP = new vector<int>[E_length];//1+
	vector<int>* E2P = new vector<int>[E_length+1];
	for (int i = 0; i < n; i++) {
		EtoP[PtoE_index[i]].emplace_back(i);
		if(PtoE_value[i] <= sqr_r){
			E2P[PtoE_index[i]].emplace_back(i);
		}
		else{
			E2P[E_length].emplace_back(i);
		}
	}
	if (E2P[E_length].size() > 0) {
		outliers = true;
	}
	//cout<< E2P[E_length].size() << endl;

	delete[] PtoE_index;
	delete[] PtoE_value;
	//PtoE_value = NULL;

	//tt = clock();

	//construct Rtree for each Core
	int npsize = 0;
	vector<Point* > near_point;
	Rtree *rtree = new Rtree[E_length+1];//1+
	for (int i = 0; i < E_length+1; ++i) {//+1
		npsize = E2P[i].size();
		near_point.resize(npsize);
		for (int j = 0; j < npsize; ++j ) {
			near_point[j] = points[E2P[i][j]];
		}
		rtree[i].bulkLoadRtree(near_point);

		vector<Point*>().swap(near_point);
		near_point.clear();
	}


	//double rq = 0;


	//label Core points
	//int factor = 2;
	//double sqr_bound = (factor * r + Eps) * (factor * r + Eps);
	double sqr_bound = (2 * r + Eps) * (2 * r + Eps);
	double sqr_eps = Eps * Eps;

	Point* cur;
	int ccount = 0;
	int epsize = 0;
	vector<Point*> targetPlace;
	vector<Point*> SumField;
	vector<int> E_neighbor;
	//start loop
	for (int i = 0; i < E_length; ++i) {
		//store Near Core index
		for (int j = 0; j < E_length; ++j) {
			if (query_distE(dist_E, i, j, E_length) <= sqr_bound) {
				E_neighbor.emplace_back(j);
			}
		}
		if (outliers) {
			E_neighbor.emplace_back(E_length);
		}
		//E_neighbor.emplace_back(E_length);

		//Foreach in this Core
		epsize = EtoP[i].size();
		for (int t = 0; t < epsize; ++t) {
			ccount = 0;
			cur = points[EtoP[i][t]];


			//t3 = clock();

			for(vector<int>::iterator iter = E_neighbor.begin() ; iter!=E_neighbor.end() ; ++iter){
				rtree[*iter].rangeQuery(cur,sqr_eps,targetPlace);
				ccount += targetPlace.size();
				SumField.reserve(ccount);
				SumField.insert(SumField.end(),targetPlace.begin(),targetPlace.end());
				vector<Point*>().swap(targetPlace);
				targetPlace.clear();
			}




			//when it is corepoint
			if (ccount >= MinPts) {
				cur->isCorePoint = true;
				cur->countn = ccount;
				//points[i]->neighbor = &this_neighbor[0];
				cur->neighbor = new Point * [ccount];
				for (int k = 0; k < ccount; k++) {
					//points[index]->neighbor[t] = this_neighbor[t];
					cur->neighbor[k] = SumField[k];
				}
			}

			vector<Point*>().swap(SumField);
			SumField.clear();
		}

		vector<int>().swap(E_neighbor);
		E_neighbor.clear();
	}



	for (int i = 0; i < E_length+1; ++i) {
		rtree[i].releaseSpace();
	}


	//cout <<"rangequery time\t"<< rq / (double)CLOCKS_PER_SEC << endl;



	//release memory
	delete[] dist_E;
	dist_E = NULL;
	for (int i = 0; i < E_length; i++) {
		//EtoP[i].clear();
		vector<int>().swap(EtoP[i]);
		vector<int>().swap(E2P[i]);
	}
	vector<int>().swap(E2P[E_length]);
	delete[] EtoP;
	delete[] E2P;
	EtoP = NULL;
	E2P = NULL;

}*/

/*void CorePoint(vector<Point*>& points, int* PtoE_index, unsigned int* PtoE_value, vector<int>& E, int E_length, int MinPts, double Eps, double r, int n, int d) {

	bool outliers = false;

	//calculate dist_E
	unsigned int* dist_E = new unsigned int[(E_length - 1) * E_length / 2];
	generate_distE(dist_E, points, E, E_length, d);


	double sqr_r = r * r;
	//calculate EtoP
	//vector<int>* EtoP = new vector<int>[E_length];//1+
	vector<int>* E2P = new vector<int>[E_length + 1];
	for (int i = 0; i < n; i++) {
		//EtoP[PtoE_index[i]].emplace_back(i);
		if (PtoE_value[i] <= sqr_r) {
			E2P[PtoE_index[i]].emplace_back(i);
		}
		else {
			E2P[E_length].emplace_back(i);
		}
	}
	if (E2P[E_length].size() > 0) {
		outliers = true;
	}
	//cout<< E2P[E_length].size() << endl;

	//delete[] PtoE_index;
	delete[] PtoE_value;
	//PtoE_value = NULL;

	//tt = clock();

	//construct Rtree for each Core
	int npsize = 0;
	vector<Point* > near_point;
	Rtree* rtree = new Rtree[E_length + 1];//1+
	for (int i = 0; i < E_length + 1; ++i) {//+1
		npsize = E2P[i].size();
		near_point.resize(npsize);
		for (int j = 0; j < npsize; ++j) {
			near_point[j] = points[E2P[i][j]];
		}
		rtree[i].bulkLoadRtree(near_point);

		vector<Point*>().swap(near_point);
		near_point.clear();
	}


	//double rq = 0;


	//label Core points
	//int factor = 2;
	//double sqr_bound = (factor * r + Eps) * (factor * r + Eps);
	double sqr_bound = (2 * r + Eps) * (2 * r + Eps);
	double sqr_eps = Eps * Eps;

	Point* cur;
	int ccount = 0;
	int epsize = 0;
	vector<Point*> targetPlace;
	vector<Point*> SumField;
	vector<int> E_neighbor;
	//start loop
	for (int i = 0; i < E_length; ++i) {
		//store Near Core index
		for (int j = 0; j < E_length; ++j) {
			if (query_distE(dist_E, i, j, E_length) <= sqr_bound) {
				E_neighbor.emplace_back(j);
			}
		}
		if (outliers) {
			E_neighbor.emplace_back(E_length);
		}
		//E_neighbor.emplace_back(E_length);

		//Foreach in this Core
		epsize = E2P[i].size();
		for (int t = 0; t < epsize; ++t) {
			ccount = 0;
			cur = points[E2P[i][t]];


			//t3 = clock();

			for (vector<int>::iterator iter = E_neighbor.begin(); iter != E_neighbor.end(); ++iter) {
				rtree[*iter].rangeQuery(cur, sqr_eps, targetPlace);
				ccount += targetPlace.size();
				SumField.reserve(ccount);
				SumField.insert(SumField.end(), targetPlace.begin(), targetPlace.end());
				vector<Point*>().swap(targetPlace);
				targetPlace.clear();
			}




			//when it is corepoint
			if (ccount >= MinPts) {
				cur->isCorePoint = true;
				cur->countn = ccount;
				//points[i]->neighbor = &this_neighbor[0];
				cur->neighbor = new Point * [ccount];
				for (int k = 0; k < ccount; k++) {
					//points[index]->neighbor[t] = this_neighbor[t];
					cur->neighbor[k] = SumField[k];
				}
			}

			vector<Point*>().swap(SumField);
			SumField.clear();
		}


		vector<int>().swap(E_neighbor);
		E_neighbor.clear();
	}


	if (outliers) {
		double sqr_b = (r + Eps) * (r + Eps);
		int epsize = E2P[E_length].size();
		int curind;
		for (int i = 0; i < epsize; i++) {
			ccount = 0;
			curind = E2P[E_length][i];
			cur = points[curind];

			for (int j = 0; j < E_length; ++j) {
				if (square_dist(points, curind, E[j], d) <= sqr_b) {
					E_neighbor.emplace_back(j);
				}
			}
			E_neighbor.emplace_back(E_length);

			for (vector<int>::iterator iter = E_neighbor.begin(); iter != E_neighbor.end(); ++iter) {
				rtree[*iter].rangeQuery(cur, sqr_eps, targetPlace);
				ccount += targetPlace.size();
				SumField.reserve(ccount);
				SumField.insert(SumField.end(), targetPlace.begin(), targetPlace.end());
				vector<Point*>().swap(targetPlace);
				targetPlace.clear();
			}




			//when it is corepoint
			if (ccount >= MinPts) {
				cur->isCorePoint = true;
				cur->countn = ccount;
				//points[i]->neighbor = &this_neighbor[0];
				cur->neighbor = new Point * [ccount];
				for (int k = 0; k < ccount; k++) {
					//points[index]->neighbor[t] = this_neighbor[t];
					cur->neighbor[k] = SumField[k];
				}
			}

			vector<Point*>().swap(SumField);
			SumField.clear();
			vector<int>().swap(E_neighbor);
			E_neighbor.clear();
		}
	}


	for (int i = 0; i < E_length + 1; ++i) {
		rtree[i].releaseSpace();
	}


	//cout <<"rangequery time\t"<< rq / (double)CLOCKS_PER_SEC << endl;



	//release memory
	delete[] dist_E;
	dist_E = NULL;
	for (int i = 0; i < E_length; i++) {
		//EtoP[i].clear();
		//vector<int>().swap(EtoP[i]);
		vector<int>().swap(E2P[i]);
	}
	vector<int>().swap(E2P[E_length]);
	//delete[] EtoP;
	delete[] E2P;
	//EtoP = NULL;
	E2P = NULL;
	delete[] PtoE_index;
}*/

void CorePoint(vector<Point *> &points, int *PtoE_index, unsigned int *PtoE_value, vector<int> &E, int E_length, int MinPts, double Eps, double r, int n, int d)
{

	clock_t t1 = clock();
	bool outliers = false;

	// calculate dist_E
	unsigned int *dist_E = new unsigned int[(E_length - 1) * E_length / 2];
	generate_distE(dist_E, points, E, E_length, d);

	double sqr_r = r * r;
	// calculate EtoP
	// vector<int>* EtoP = new vector<int>[E_length];//1+
	vector<int> *E2P = new vector<int>[E_length + 1];
	for (int i = 0; i < n; i++)
	{
		// EtoP[PtoE_index[i]].emplace_back(i);
		if (PtoE_value[i] <= sqr_r)
		{
			E2P[PtoE_index[i]].emplace_back(i);
		}
		else
		{
			E2P[E_length].emplace_back(i);
		}
	}
	if (E2P[E_length].size() > 0)
	{
		outliers = true;
	}
	// cout<< E2P[E_length].size() << endl;

	// delete[] PtoE_index;
	delete[] PtoE_value;
	// PtoE_value = NULL;

	clock_t t2 = clock();
	cout << "t2-t1:\t" << (t2 - t1) / (double)CLOCKS_PER_SEC << endl;

	// construct Rtree for each Core
	int npsize = 0;
	vector<Point *> near_point;
	Rtree *rtree = new Rtree[E_length + 1]; // 1+
	for (int i = 0; i < E_length + 1; ++i)
	{ //+1
		npsize = E2P[i].size();
		near_point.resize(npsize);
		for (int j = 0; j < npsize; ++j)
		{
			near_point[j] = points[E2P[i][j]];
		}
		rtree[i].bulkLoadRtree(near_point);

		vector<Point *>().swap(near_point);
		near_point.clear();
	}

	// double rq = 0;

	// label Core points
	// int factor = 2;
	// double sqr_bound = (factor * r + Eps) * (factor * r + Eps);
	double sqr_bound = (2 * r + Eps) * (2 * r + Eps);
	double sqr_eps = Eps * Eps;

	Point *cur;
	int ccount = 0;
	int epsize = 0;
	vector<Point *> targetPlace;
	vector<Point *> SumField;
	vector<int> E_neighbor;
	// start loop
	for (int i = 0; i < E_length; ++i)
	{
		// store Near Core index
		for (int j = 0; j < E_length; ++j)
		{
			if (query_distE(dist_E, i, j, E_length) <= sqr_bound)
			{
				E_neighbor.emplace_back(j);
			}
		}
		if (outliers)
		{
			E_neighbor.emplace_back(E_length);
		}
		// E_neighbor.emplace_back(E_length);

		// Foreach in this Core
		epsize = E2P[i].size();
		for (int t = 0; t < epsize; ++t)
		{
			ccount = 0;
			cur = points[E2P[i][t]];

			// t3 = clock();

			for (vector<int>::iterator iter = E_neighbor.begin(); iter != E_neighbor.end(); ++iter)
			{
				rtree[*iter].rangeQuery(cur, sqr_eps, targetPlace);
				ccount += targetPlace.size();
				SumField.reserve(ccount);
				SumField.insert(SumField.end(), targetPlace.begin(), targetPlace.end());
				vector<Point *>().swap(targetPlace);
				targetPlace.clear();
			}

			// when it is corepoint
			if (ccount >= MinPts)
			{
				cur->isCorePoint = true;
				cur->countn = ccount;
				// points[i]->neighbor = &this_neighbor[0];
				cur->neighbor = new Point *[ccount];
				for (int k = 0; k < ccount; k++)
				{
					// points[index]->neighbor[t] = this_neighbor[t];
					cur->neighbor[k] = SumField[k];
				}
			}

			vector<Point *>().swap(SumField);
			SumField.clear();
		}

		vector<int>().swap(E_neighbor);
		E_neighbor.clear();
	}

	if (outliers)
	{
		double sqr_b = (r + Eps) * (r + Eps);
		double sqr_bd = 4 * (r + Eps) * (r + Eps);
		int epsize = E2P[E_length].size();
		int curind;
		int eind;
		for (int i = 0; i < epsize; i++)
		{
			ccount = 0;
			curind = E2P[E_length][i];
			cur = points[curind];
			eind = PtoE_index[curind];
			for (int j = 0; j < E_length; ++j)
			{
				if (query_distE(dist_E, eind, j, E_length) <= sqr_bd)
				{
					if (square_dist(points, curind, E[j], d) <= sqr_b)
					{
						E_neighbor.emplace_back(j);
					}
				}
			}
			E_neighbor.emplace_back(E_length);

			for (vector<int>::iterator iter = E_neighbor.begin(); iter != E_neighbor.end(); ++iter)
			{
				rtree[*iter].rangeQuery(cur, sqr_eps, targetPlace);
				ccount += targetPlace.size();
				SumField.reserve(ccount);
				SumField.insert(SumField.end(), targetPlace.begin(), targetPlace.end());
				vector<Point *>().swap(targetPlace);
				targetPlace.clear();
			}

			// when it is corepoint
			if (ccount >= MinPts)
			{
				cur->isCorePoint = true;
				cur->countn = ccount;
				// points[i]->neighbor = &this_neighbor[0];
				cur->neighbor = new Point *[ccount];
				for (int k = 0; k < ccount; k++)
				{
					// points[index]->neighbor[t] = this_neighbor[t];
					cur->neighbor[k] = SumField[k];
				}
			}

			vector<Point *>().swap(SumField);
			SumField.clear();
			vector<int>().swap(E_neighbor);
			E_neighbor.clear();
		}
	}

	for (int i = 0; i < E_length + 1; ++i)
	{
		rtree[i].releaseSpace();
	}

	// cout <<"rangequery time\t"<< rq / (double)CLOCKS_PER_SEC << endl;

	// release memory
	delete[] dist_E;
	dist_E = NULL;
	for (int i = 0; i < E_length; i++)
	{
		// EtoP[i].clear();
		// vector<int>().swap(EtoP[i]);
		vector<int>().swap(E2P[i]);
	}
	vector<int>().swap(E2P[E_length]);
	// delete[] EtoP;
	delete[] E2P;
	// EtoP = NULL;
	E2P = NULL;
	delete[] PtoE_index;
}

// 聚类
void DBSCAN(vector<Point *> &points, int n)
{
	// bool* visited = new bool[n];
	vector<Point *> local;
	int newlabel = 1;
	// int res;
	// int it;
	// memset(visited, false, n);

	// bool* added = new bool[n];
	// memset(added, false, n);

	for (int i = 0; i < n; ++i)
	{
		if (points[i]->visited)
		{
			continue;
		}
		if (points[i]->isCorePoint)
		{
			points[i]->clustringID = newlabel;
			points[i]->visited = true;
			local.emplace_back(points[i]);
			// added[i] = true;

			for (int j = 0; j < local.size(); ++j)
			{
				// res = local[j];

				for (int t = 0; t < local[j]->countn; ++t)
				{
					// local[j]->neighbor[t];

					if (!local[j]->neighbor[t]->visited)
					{
						local[j]->neighbor[t]->visited = true;
						local[j]->neighbor[t]->clustringID = newlabel;

						// if (points[it]->isCorePoint && !added[it]) {
						if (local[j]->neighbor[t]->isCorePoint)
						{
							local.emplace_back(local[j]->neighbor[t]);
							// added[it] = true;
						}
					}
				}
			}
			newlabel++;
			vector<Point *>().swap(local);
			local.clear();
		}
	}
	cata = newlabel - 1;
	// delete[] visited;
	// visited = NULL;

	// delete[] added;
	// added = NULL;
}

// 初始化PtoE
void init_PtoE(vector<Point *> &pt, int *PtoE_index, unsigned int *PtoE_value, int *E_new, int E_length, int act_sel, int n, int d)
{
	unsigned int min;
	unsigned int current;
	for (int i = 0; i < n; ++i)
	{
		if (pt[i]->fromE)
		{
			continue;
		}
		min = square_dist(pt, i, E_new[0], d);
		PtoE_index[i] = E_length;
		for (int j = 1; j < act_sel; ++j)
		{
			current = square_dist(pt, i, E_new[j], d);
			if (current < min)
			{
				min = current;
				PtoE_index[i] = E_length + j;
			}
		}
		PtoE_value[i] = min;
	}
}

// 更新PtoE
void update_PtoE(vector<Point *> &pt, int *PtoE_index, unsigned int *PtoE_value, int *E_new, int E_length, int act_sel, int n, int d)
{
	unsigned int min;
	unsigned int current;
	for (int i = 0; i < n; ++i)
	{
		if (pt[i]->fromE)
		{
			continue;
		}
		min = PtoE_value[i];
		for (int j = 0; j < act_sel; ++j)
		{
			current = square_dist(pt, i, E_new[j], d);
			if (current < min)
			{
				min = current;
				PtoE_index[i] = E_length + j;
			}
		}
		PtoE_value[i] = min;
	}
}

// add to E, remove from P.
void AddtoE(vector<Point *> &pt, vector<int> &E, int *PtoE_index, unsigned int *PtoE_value, int *E_new, int E_length, int act_sel)
{
	E.resize(E.size() + act_sel); // cout << act_sel << endl;
	for (int i = 0; i < act_sel; ++i)
	{
		// 移除P
		pt[E_new[i]]->fromE = true;
		// 加入E
		E[E_length + i] = E_new[i];
		// save EtoE
		PtoE_index[E_new[i]] = E_length + i;
		PtoE_value[E_new[i]] = 0;
	}
	// cout << act_sel << endl;//E.assign(E_new,E_new+act_sel+1);
}

unsigned int square_dist(vector<Point *> &pt, int a, int b, int d)
{
	// if (a > n || b > n) {
	//	cout << "" << endl;
	// }

	register unsigned int temp;
	// double temp;
	register unsigned int dist = 0; // sue long/longlong, when distances are large
	// double dist = 0;
	for (int i = 0; i < d; ++i)
	{
		temp = pt[a]->coords[i] - pt[b]->coords[i];
		dist += temp * temp;
	}
	return dist;
}

// 调整堆
void AdjustDown(int *top, unsigned int *topvalue, int root, int size)
{
	int child = 2 * root + 1; // 下标从0开始,左孩子
	while (root < size / 2)
	{
		if (child + 1 < size && topvalue[child] > topvalue[child + 1])
		{
			child++;
		}
		if (child < size && topvalue[root] > topvalue[child])
		{
			swap(topvalue[root], topvalue[child]);
			swap(top[root], top[child]);
			root = child;
			child = 2 * root + 1;
		}
		else
		{
			break;
		}
	}
}

// 计算最远k个点
unsigned int furthest_KinP(int *Q, vector<Point *> &pt, unsigned int *PtoE_value, int k, int n)
{

	// 堆排序法
	unsigned int *Qvalue = new unsigned int[k];
	int j = 0;
	for (int i = 0; i < k; ++i)
	{
		while (j < n && pt[j]->fromE)
		{
			j++;
		}
		if (j >= n)
		{
			cout << "error furthest_KinP\n";
			k = i;
			exit(0);
			break;
		}
		Qvalue[i] = PtoE_value[j];
		Q[i] = j;
		j++;
	}
	for (int i = k / 2 - 1; i >= 0; i--)
	{
		AdjustDown(Q, Qvalue, i, k);
	}
	for (; j < n; j++)
	{
		if (PtoE_value[j] > Qvalue[0])
		{
			Qvalue[0] = PtoE_value[j];
			Q[0] = j;
			AdjustDown(Q, Qvalue, 0, k);
		}
	}
	unsigned int dis_QtoE = Qvalue[0];
	delete[] Qvalue;
	Qvalue = NULL;
	return dis_QtoE;
}

// 随机选择0~sup-1中的num个不同数字
void select_random_1(int *pick, int num, int sup)
{
	if (sup < num)
	{
		cout << "error select_random" << endl;
		exit(0);
		return;
	}
	int temp;
	int flag;
	srand((unsigned)(time(0)));
	// 总数很大适合选用
	for (int i = 0; i < num; ++i)
	{
		while (true)
		{
			flag = 0;
			temp = (rand() % sup);
			for (int j = 0; j < i; ++j)
			{
				if (pick[j] == temp)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
			{
				pick[i] = temp;
				break;
			}
		}
	}
}

void select_random_2(int *pick, int num, int sup)
{
	if (sup < num)
	{
		cout << "error select_random" << endl;
		return;
	}
	int temp;
	srand((unsigned)(time(0)));

	int *line = new int[sup];
	for (int i = 0; i < sup; ++i)
	{
		line[i] = i;
	}
	for (int i = 0; i < num; ++i)
	{
		temp = (rand() % (sup - i)) + i;
		pick[i] = line[temp];
		line[temp] = line[i];
	}
	delete[] line;
	line = NULL;
}

// 删除E_new相近的点
int delete_nearE(vector<Point *> &pt, int *E_new, int beforDel, double r, int d)
{
	int afterDel = beforDel;
	for (int i = 1; i < afterDel; ++i)
	{
		for (int j = 0; j < i; j++)
		{
			if (square_dist(pt, E_new[i], E_new[j], d) < r)
			{
				E_new[i] = E_new[afterDel - 1];
				afterDel--;
				i--;
				break;
			}
		}
	}
	return afterDel;
}

// 生成E间距离
void generate_distE(unsigned int *dist_E, vector<Point *> &points, vector<int> &E, int E_length, int d)
{
	// 三角矩阵
	int index;
	for (int i = 0; i < E_length - 1; ++i)
	{
		for (int j = i + 1; j < E_length; ++j)
		{
			index = i * E_length - (3 + i) * i / 2 + j - 1;
			dist_E[index] = square_dist(points, E[i], E[j], d);
		}
	}
	// 二维数组
}

// 查询E间距离
inline unsigned int query_distE(unsigned int *dist_E, int a, int b, int E_length)
{
	if (a == b)
	{
		return 0;
	}
	int index;
	if (a > b)
	{
		// swap(a, b);
		index = b * E_length - (3 + b) * b / 2 + a - 1;
	}
	else
	{
		index = a * E_length - (3 + a) * a / 2 + b - 1;
	}
	return dist_E[index];
}

// 读取数据
void readDatasetFromFile(char *filePath, vector<Point *> &pt, int n, int d)
{

	FILE *file = fopen(filePath, "rt");
	pt.clear();
	// pt.reserve(n);
	pt.resize(n);
	for (int i = 0; i < n; ++i)
	{
		pt[i] = readPointFromFile(file, d);
	}
	fclose(file);
}

Point *readPointFromFile(FILE *fileHandle, int d)
{

	int *coords = new int[d];
	for (int i = 0; i < d; ++i)
	{
		fscanf(fileHandle, "%d", &(coords[i]));
	}
	Point *pt = new Point(coords);
	return pt;
}

// 写入文件
void writeClusterToFile(vector<Point *> &points, char *filePath, int n)
{
	ofstream myout(filePath);
	for (int i = 0; i < n; ++i)
	{
		myout << points[i]->clustringID << endl;
	}
	myout.close();
}
