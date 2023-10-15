#include "../header/DBSCAN.h"
#include <getopt.h>
using namespace std;
int main(int argc, char *argv[])
{
    cout << argc << " args received." << endl;
    int o;
    int alg;
    int dim;          // dimension
    int n;            // dataset scale
    int MinPts;       // MinPts
    double epsilon;   // epsilon
    double rho = 0.5; // ρ

    int init_num;
    int farthest_num;
    int select_num;

    int max_layer;
    int terminate_layer;

    string input_fileName;
    string output_data_filename;
    string output_center_filename;

    static const struct option longopts[] = {
        {"alg", required_argument, NULL, 'a'},
        {"dim", required_argument, NULL, 'd'},
        {"n", required_argument, NULL, 'n'},
        {"minpts", required_argument, NULL, 'm'},
        {"epsilon", required_argument, NULL, 'e'},
        {"rho", required_argument, NULL, 'r'},
        {"init_num", required_argument, NULL, 'i'},
        {"farthest_num", required_argument, NULL, 'f'},
        {"select_num", required_argument, NULL, 's'},
        {"max_layer", required_argument, NULL, 'u'},
        {"terminate_layer", required_argument, NULL, 'v'},
        {"input", required_argument, NULL, 'x'},
        {"output_data", required_argument, NULL, 'y'},
        {"output_center", required_argument, NULL, 'z'},
        {NULL, 0, NULL, 0}};
    while ((o = getopt_long(argc, argv, "a:d:n:m:e:r:i:f:s:u:v:x:y:z:", longopts, NULL)) != -1)
    {
        switch (o)
        {
        case 'a':
            alg = atoi(optarg);
            break;
        case 'd':
            dim = atoi(optarg);
            break;
        case 'n':
            n = atoi(optarg);
            break;
        case 'm':
            MinPts = atoi(optarg);
            break;
        case 'e':
            epsilon = atof(optarg);
            break;
        case 'r':
            rho = atof(optarg);
            break;
        case 'i':
            init_num = atoi(optarg);
            break;
        case 'f':
            farthest_num = atoi(optarg);
            break;
        case 's':
            select_num = atoi(optarg);
            break;
        case 'u':
            max_layer = atoi(optarg);
            break;
        case 'v':
            terminate_layer = atoi(optarg);
            break;
        case 'x':
            input_fileName = optarg;
            break;
        case 'y':
            output_data_filename = optarg;
            break;
        case 'z':
            output_center_filename = optarg;
            break;
        }
    }
    cout << "Running algorithm: " << alg << endl;
    switch (alg)
    {
    // alg1:vanilla DBSCAN
    // alg2:our exact DBSCAN
    // alg3:our exact DBSCAN via CoverTree
    // alg4:our exact DBSCAN with edit distance
    // alg5:our approximate DBSCAN
    // alg6:our approximate DBSCAN via Randomized Gonzalez
    // alg7:our approximate DBSCAN with edit distance
    case 1:
        DBSCAN_vanilla(dim, n, MinPts, epsilon, input_fileName, output_data_filename, output_center_filename);
        break;
    case 2:
        DBSCAN_exact(dim, n, MinPts, epsilon, rho, input_fileName, output_data_filename, output_center_filename);
        break;
    case 3:
        DBSCAN_exact_covertree(dim, n, MinPts, epsilon, rho, input_fileName, output_data_filename, output_center_filename, max_layer, terminate_layer);
        break;
    case 4:
        DBSCAN_exact_edit(dim, n, MinPts, epsilon, rho, input_fileName, output_data_filename, output_center_filename);
        break;
    case 5:
        DBSCAN_approx(dim, n, MinPts, epsilon, rho, input_fileName, output_data_filename, output_center_filename);
        break;
    case 6:
        DBSCAN_approx_random(dim, n, MinPts, epsilon, rho, init_num, farthest_num, select_num, input_fileName, output_data_filename, output_center_filename);
        break;
    case 7:
        DBSCAN_approx_edit(dim, n, MinPts, epsilon, rho, input_fileName, output_data_filename, output_center_filename);
        break;
    case 8:
        DBSCAN_approx_covertree(dim, n, MinPts, epsilon, rho, input_fileName, output_data_filename, output_center_filename, max_layer, terminate_layer);
        break;
    case 9:
        DBSCAN_vanilla_edit(dim, n, MinPts, epsilon, input_fileName, output_data_filename, output_center_filename);
        break;
        /* case 10:
        DBSCAN_pp_random(dim, n, MinPts, epsilon, input_fileName, output_data_filename, output_center_filename) case 11:
            DBSCAN_pp_kcenter(dim, n, MinPts, epsilon, input_fileName, output_data_filename, output_center_filename) */
    }

    int terminate = 1;
    exit(0);
    return 0;
}