#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "multinomial.h"
#include <boost/math/special_functions/binomial.hpp>

using namespace std;

double Prob(double L, uint n, double p, uint k)
{

    vector<double> readPr;
    for (int i = 0; i <= k; ++i)
        readPr.push_back(boost::math::binomial_coefficient<long double>(k, i) * pow(double(n) / L, i) * pow(double(L - n) / L, k - i));

//    for (auto s : readPr)
//        cout << s << endl;
    // Вектор, в котором i-ое значение  - это вероятность ошибки при попадании в конкретный
    // нуклиотид i раз
    vector<double> FaultPr;
    FaultPr.push_back(0.75);
    for (int j = 1; j <= k; ++j)
    {
        if (p != 0)
        {
            double fault_tmp = 0;
            for (long unsigned int i = 0; i <= j; ++i)
            {
                for (long unsigned int s = 0; s <= j; ++s)
                {
                    for (long unsigned int r = 0; r <= j; ++r)
                    {
                        for (long unsigned int e = 0; e <= j; ++e)
                        {
                            if (i + s + r + e == j)
                            {
                                if ((i > s && i > r && i > e) || (s > i && s > r && s > e) || (r > i && r > s && r > e))
                                    fault_tmp += multinomial::multi<long double>({i, s, r, e}) * pow(p / 3, i) * pow(p / 3, s) * pow(p / 3, r) * pow(1 - p, e);
                                else if ((i == e && i > r && i > s) || (s == e && s > i && s > r) || (r == e && r > s && r > i))
                                    fault_tmp += multinomial::multi<long double>({i, s, r, e}) * pow(p / 3, i) * pow(p / 3, s) * pow(p / 3, r) * pow(1 - p, e) * 0.5;
                                else if ((i == e && s == e && i > r) || (s == e && r == e && s > i) || (r == e && i == e && r > s))
                                    fault_tmp += multinomial::multi<long double>({i, s, r, e}) * pow(p / 3, i) * pow(p / 3, s) * pow(p / 3, r) * pow(1 - p, e) * 2 / 3.0;
                                else if (i == e && s == e && r == e)
                                    fault_tmp += multinomial::multi<long double>({i, s, r, e}) * pow(p / 3, i) * pow(p / 3, s) * pow(p / 3, r) * pow(1 - p, e) * 0.75;
                            }
                        }
                    }
                }
            }
            FaultPr.push_back(fault_tmp);
        }
        else
        {
            FaultPr.push_back(0.0);
        }
    }

    double result = 0;
    for (int i = 0; i <= k; ++i)
        result += readPr[i] * FaultPr[i];

    return result * L;
}

int main()
{
    double L;
    uint n;
    double p;
    uint k;

    string line;
    ifstream tests_file;
    ofstream output;
    tests_file.open("2.txt", ios_base::in);
    output.open("output.txt");
    cout.precision(10);
    output.precision(10);

    while (getline(tests_file, line))
    {
        istringstream ss(line);
        ss >> L >> n >> p >> k;
        output << fixed << Prob(L, n, p, k) << endl;
    };

    output.close();
    tests_file.close();
}