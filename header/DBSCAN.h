#pragma once
#include <queue>
#include <map>
#include "Kcenter.h"
#include "CoverTree.h"
void range_search(vector<Point *> &data, int MinPts, double epsilon, vector<vector<double>> distance_matrix);
void range_search_edit(vector<Point *> &data, int MinPts, double epsilon, vector<vector<double>> distance_matrix);
int cluster_core_vanilla(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix);
void calculate_epsilon_net(CoverTree *T, int epsilon_net_layer, vector<vector<double>> &distance_matrix, vector<Center *> &centers);
void calculate_Ap(vector<Center *> &centers, double epsilon, vector<vector<double>> &distance_matrix);
void calculate_Ap_covertree(vector<Center *> &centers, vector<Point *> &data, double r_max, double epsilon, vector<vector<double>> &distance_matrix);
void mark_core(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C);
void mark_core_edit(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C);
void mark_core_covertree(vector<Center *> &centers, CoverTree *T, vector<vector<double>> &distance_matrix, int minpts, double epsilon, vector<Point *> &data, vector<Center *> &C);
Node *NN_search(Node *p, Node *root, vector<vector<double>> &distance_matrix, vector<Point *> &data, int terminate_layer);
Node *NN_search2(Node *p, Center *root, vector<vector<double>> &distance_matrix, vector<Point *> &data, int terminate_layer);
void connect_centers(vector<Center *> &centers, vector<vector<double>> &distance_matrix, vector<Point *> &data, int terminate_layer, double r_max, int minpts, double epsilon);
void connect_centers_bruteforce(vector<Center *> &centers, vector<vector<double>> &distance_matrix, vector<Point *> &data, int terminate_layer, double r_max, int minpts, double epsilon);
double BCP(vector<Point *> &c1, vector<Point *> &c2, double epsilon, vector<vector<double>> &distance_matrix);
double BCP_edit(vector<Point *> &c1, vector<Point *> &c2, double epsilon, vector<vector<double>> &distance_matrix);
int cluster_core(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C);
int cluster_core_approx(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C);
int cluster_core_edit(vector<Point *> &data, vector<Center *> &centers, double epsilon, int MinPts, vector<vector<double>> &distance_matrix, vector<Center *> &C);
void cluster_border(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix, vector<Center *> &C);
void cluster_border_edit(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix, vector<Center *> &C);
void cluster_border_approx(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix, vector<Center *> &C);
void cluster_border_approx_edit(vector<Point *> &data, double epsilon, vector<vector<double>> &distance_matrix, vector<Center *> &C);
void DBSCAN_vanilla(int dim, int n, int MinPts, double epsilon, string input_fileName, string output_data_filename, string output_center_filename);
void DBSCAN_vanilla_edit(int dim, int n, int MinPts, double epsilon, string input_fileName, string output_data_filename, string output_center_filename);
void DBSCAN_exact(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename);
void DBSCAN_exact_covertree(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename, int max_layer, int terminate_layer);
void DBSCAN_exact_edit(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename);
void DBSCAN_approx(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename);
void DBSCAN_approx_random(int dim, int n, int MinPts, double epsilon, double rho, int init_num, int farthest_num, int select_num, string input_fileName, string output_data_filename, string output_center_filename);
void DBSCAN_approx_edit(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename);
void DBSCAN_approx_covertree(int dim, int n, int MinPts, double epsilon, double rho, string input_fileName, string output_data_filename, string output_center_filename, int max_layer, int terminate_layer);
