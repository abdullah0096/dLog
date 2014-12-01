/* 
 * File:   main.cpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:24 AM
 */

#include "discreteLog.hpp"
#include <NTL/vec_ZZ.h>

using namespace std;

int main(int argc, char** argv) {

    long double time = 0.0;
    int end = 1;

    for (long int i = 0; i < end; i++) {

        long long int n, r, l, t;
        ZZ p, orderOfG;
        ZZ_pX g, h, irrdPoly;

        ifstream fin("in.txt");
        if (!fin) {
            cout << "\n ERROR in Main reading File in.txt...\n";
            exit(1);
        }

        cout << "\n Reading Input from in.txt\n";
        fin >> p >> n >> r >> orderOfG >> l >>t;
        ZZ_p::init(p);

        fin >> g;
        cout << "\n g:: " << g << endl;
        fin>> h;
        cout << "\n h:: " << h << endl;
        fin>>irrdPoly;
        cout << "\n irrdPoly:: " << irrdPoly << endl;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\n\n";

        discreteLog DLP(p, n, r, l, g, h, irrdPoly, t, orderOfG);
        DLP.printParameters();
        DLP.printMultipliers();

        DLP.bruteForceDL();
        //        DLP.teskeDL();

        if (DLP.cheonDL4() == 0) {
            cout << "\n Something Went Wrong.....\n";
        } else {
            time += DLP.tableGenerationTime;
            cout << "\n in Main :: Iteration # " << i << endl;
        }

        fin.close();
    }

    cout << "\n TOTAL time :: " << time / end << endl;
    return 0;
}