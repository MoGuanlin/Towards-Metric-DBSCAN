#include "../header/IO.h"
void load_data(vector<Point *> &data, int dim, int n, string &input_fileName)
{
    // initialization
    for (int i = 0; i < n; i++)
    {
        Point *newpoint = new Point(dim);
        newpoint->ID = i;
        data.push_back(newpoint);
    }
    // read coords from file
    fstream fin;
    fin.open(input_fileName, ios::in);
    if (!fin)
    {
        cout << "File not exists!" << endl;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            double value;
            fin >> value;
            data[i]->set_coords(j, value);
        }
    }
}
void load_text_data(vector<Point *> &data, int n, string &input_fileName)
{
    for (int i = 0; i < n; i++)
    {
        Point *newpoint = new Point();
        newpoint->ID = i;
        data.push_back(newpoint);
    }
    // read coords from file
    fstream fin;
    fin.open(input_fileName, ios::in);
    if (!fin)
    {
        cout << "File not exists!" << endl;
    }
    for (int i = 0; i < n; i++)
    {
        getline(fin, data[i]->text);
    }
}
void output(vector<Point *> &data, vector<Center *> &centers, string output_data_filename, string output_center_filename)
{
    ofstream f;
    ofstream f_center;
    f.open(output_data_filename, ios::trunc);
    f_center.open(output_center_filename, ios::trunc);
    for (auto it = data.begin(); it != data.end(); it++)
    {
        // output result
        f << (*it)->clusterID << '\n';
    }
    /*     for (auto it = centers.begin(); it != centers.end(); it++)
        {
            for (auto it2 = ((*it)->coords).begin(); it2 != ((*it)->coords).end(); it2++)
            {
                f_center << *it2 << " ";
            }
            f_center << "\n";
        } */
    f.close();
    f_center.close();
}
