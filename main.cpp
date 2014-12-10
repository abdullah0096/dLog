/* 
 * File:   main.cpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:24 AM
 */
/* 
 * File:   main.cpp
 * Author: Abdullah
 *
 * Created on December 2, 2014, 3:36 PM
 */

#include <cstdlib>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <fstream>

#include "discreteLog.hpp"

using namespace NTL;
using namespace std;

int main(int argc, char** argv) {

    long double time = 0.0;
    int end = 1;

    for (long int i = 0; i < end; i++) {

        long r, l, t;
        ZZ p, n, orderOfG;
        ZZ_pX g, h, irrdPoly;

        ifstream fin("in.txt");
        if (!fin) {
            cout << "\n ERROR in Main reading File in.txt...\n";
            exit(1);
        }

        cout << "\n Reading Input from in.txt\n";
        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);

        fin >> g >> h>>irrdPoly;
        //        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
        //                << "\t l :: " << l << "\t t ::" << t << "\t";
        //        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        discreteLog DLP(p, n, r, l, g, h, irrdPoly, t, orderOfG);
        DLP.printParameters();

        DLP.cheonDL3();
        cout << "\n Time By Cheon2 ::" << DLP.getTimeByCheon() << " Seconds..." << endl;

        DLP.teske3();
        cout << "\n Time By teske ::" << DLP.getTimeByTeske() << " Seconds..." << endl;

        fin.close();
    }
    return 0;
}