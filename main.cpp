/* 
 * File:   main.cpp
 * Author: Abdullah
 *
 * Created on October 16, 2014, 11:24 AM
 */

#include "discreteLog.hpp"

using namespace std;

int main(int argc, char** argv) {

    long double time = 0.0;
    int end = 1;
    for (long int i = 0; i < end; i++) {

        long n, r, orderOfG, l, t;
        ZZ p;
        ZZ_pX g, h, irrdPoly;

        ifstream fin("in.txt");
        if (!fin) {
            cout << "\n ERROR in Main reading File in.txt...\n";
            exit(1);
        }

        cout << "\n Reading Input from in.txt\n";
        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);

        fin >> g>>h;
        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\n\n";

        discreteLog DLP(p, n, r, l, g, h, t, orderOfG);
        //        DLP.printParameters();

        if (DLP.cheonDL() == 0) {
            cout << "\n Something Went Wrong.....\n";
        } else {
            time += DLP.tableGenerationTime;
            cout << "\n in Main :: Iteration # " << i << endl;
        }
        //    cout << "\n tag :: " << DLP.getTag(g);
        ;
        //    DLP.printMultipliers();
        //    DLP.printTable();
        fin.close();
    }
    cout << "\n TOTAL time :: " << time / end << endl;
    return 0;
}