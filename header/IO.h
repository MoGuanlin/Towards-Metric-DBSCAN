#pragma once
#include "Point.h"
#include <string>
#include <fstream>
using namespace std;
void load_data(vector<Point *> &data, int dim, int n, string &input_fileName);
void output(vector<Point *> &data, vector<Center *> &centers, string output_data_filename, string output_center_filename);
void load_text_data(vector<Point *> &data, int n, string &input_fileName);