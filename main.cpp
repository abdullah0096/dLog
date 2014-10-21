/* 
 * File:   main.cpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:24 AM
 */

#include "discreteLog.hpp"

using namespace std;

int main(int argc, char** argv) {

    long long int p, n, r, orderOfG;
    ZZX g, h;

    cout << "\n Enter p, n, r, orderOfG :: ";
    cin >> p >> n >> r>>orderOfG;

    cout << "\n Enter g ,h ";
    cin >> g>>h;

    discreteLog DLP(p, n, r, 5, g, h, orderOfG);

    DLP.printParameters();
    DLP.printMultipliers();
    //    DLP.printTable();

    return 0;
}