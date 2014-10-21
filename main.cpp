/* 
 * File:   main.cpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:24 AM
 */

#include <cstdlib>
#include "discreteLog.hpp"

//#include <NTL/ZZXFactoring.h>
//using namespace NTL;

using namespace std;

int main(int argc, char** argv) {

    long long int p, n, r, orderOfG;
    //    ZZX g, h;

    cout << "\n Enter p, n, r, orderOfG :: ";
    cin >> p >> n >> r>>orderOfG;

    discreteLog DLP(p, n, r, 5, 7, 3, orderOfG);

    DLP.printParameters();
    DLP.printMultipliers();
    DLP.printTable();

    return 0;
}