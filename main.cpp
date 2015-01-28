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
#include <fstream>
#include <iomanip> 

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <NTL/GF2X.h>

#include "discreteLog.hpp"
#include "discreteLogGF2.hpp"

using namespace NTL;
using namespace std;

void foo();

int main(int argc, char** argv) {

    // <editor-fold defaultstate="collapsed" desc="GF2 Code ">
    long r, l, t;
    ZZ p, n, orderOfG;
    GF2X irrdPoly;


    long numberOfIterations = 1;

    ifstream fin("in.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    fin >>irrdPoly;
    GF2E::init(irrdPoly);
    GF2E g, h;

    cout << "\n irrdPoly :: " << irrdPoly << endl;

    cout << "\n Reading Input from in.txt\n";
    fin >> p >> n >> r >> orderOfG >> l >>t;
    ZZ_p::init(p);

    fin >> g >> h;
    cout << "\n g :: " << g << "\t h :: " << h << endl;

    long long int cnt = 0;
    GF2E tmp;

    long double time = 0;

    h = (h * h * h * h * h * h * h * h * h);
    g = (g * g * g);

    cout << "\n g * g :: " << g * g << endl;
    cout << "\n g :: " << g << endl;
    cout << "\n h :: " << h << endl;
    exit(0);
    while (1) {
        timestamp_t TableLookUpTimeStart = utility::get_timestamp();
        tmp = (h * h);
        timestamp_t TableLookUpTimeEnd = utility::get_timestamp();
        time += utility::getTimeInSeconds(TableLookUpTimeEnd, TableLookUpTimeStart);
        cnt++;
        if (cnt > 10000000)
            break;
    }
    cout << "\n time :: " << time << "\t cnt :: " << cnt << endl;
    ;
    //    discreteLogGF2 DLP(p, n, r, l, g, h, irrdPoly, t, orderOfG);
    //    DLP.printParameters();
    //    DLP.cheonDL3();
    //    cout << "\n Cheon Time :: " << DLP.getTimeByCheon() << " Seconds." << endl;
    //    cout << "\nr l  t p^n \tTime Cheon \tTable Generation Time \t Gamma Time \tInner Prod Time \tTable Look-Up Time \tMiscellaneous Time \t Actual Multiplication" << endl;
    //    cout << r << " " << l << " " << trunc(log2(r)) << " 2^" << n << "\t" << DLP.getTimeByCheon()
    //            << " Sec\t   " << DLP.getTableGenerationTime() << " Sec\t\t "
    //            << DLP.gammaTime / numberOfIterations << " Sec\t " << DLP.innerProductTime / numberOfIterations
    //            << " Sec\t\t " << DLP.tableLookUpTime / numberOfIterations << " Sec\t\t"
    //            << DLP.miscellaneousTime / numberOfIterations << " Sec\t\t " << DLP.actualMultiplicationTime / numberOfIterations << endl;
    //    cout << r << " " << l << " " << trunc(log2(r)) << " 2^" << n << "\t" << DLP.getTimeByCheon()
    //            << " Sec\t   " << DLP.getTableGenerationTime() << " Sec\t\t "
    //            << ((DLP.gammaTime / numberOfIterations) / DLP.getTimeByCheon()) * 100 << " %  \t " << ((DLP.innerProductTime / numberOfIterations) / DLP.getTimeByCheon()) * 100
    //            << " %  \t\t " << ((DLP.tableLookUpTime / numberOfIterations) / DLP.getTimeByCheon()) * 100 << " %   \t\t"
    //            << ((DLP.miscellaneousTime / numberOfIterations) / DLP.getTimeByCheon()) * 100 << " %   \t\t " << ((DLP.actualMultiplicationTime / numberOfIterations) / DLP.getTimeByCheon()) * 100 << " %" << endl;
    //
    //    cout << "\n Collision Time :: " << DLP.collisionTime << endl;
    //    cout << "\n Walk Cnt Time :: " << DLP.walkCntTime << endl;
    ;
    //    DLP.teske();
    //    cout << "\n Teske Time :: " << fixed << DLP.getTimeByTeske() << " Seconds." << endl;

    // </editor-fold>
    return 0;
}

void foo() {
    long double time = 0.0;
    long numberOfIterations = 10;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("cheon_time_test1.txt");
    ofstream teske("teske_time_test1.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        ZZ_pX g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        ZZ_p::init(p);
        fin >> g >> h>>irrdPoly;
        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;


        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLog *DLP;
        DLP = new discreteLog(p, n, r, l, g, h, irrdPoly, t, orderOfG);

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->cheonDL3();
            cheonTime += DLP->getTimeByCheon();
        }

        cheonTime = cheonTime / numberOfIterations;
        if (whileLoopCnt == 0) {
            cheon << "Number of Iterations :: " << numberOfIterations << endl;
            cheon << "Number of Iterations of Walk :: " << DLP->getNumberOfIterations() << endl;
            cheon << "\nr l  t p^n \tTime Cheon \tTable Generation Time \t Gamma Time \tInner Prod Time \tTable Look-Up Time \tMiscellaneous Time \t Actual Multiplication" << endl;
        }
        //        cheon << std::setprecision(5);

        cheon << std::fixed;
        cheon << r << " " << l << " " << trunc(log2(r)) << " 2^" << n << "\t" << cheonTime << " Sec\t   " << DLP->getTableGenerationTime() << " Sec\t\t "
                << DLP->gammaTime / numberOfIterations << " Sec\t " << DLP->innerProductTime / numberOfIterations
                << " Sec\t\t " << DLP->tableLookUpTime / numberOfIterations << " Sec\t\t"
                << DLP->miscellaneousTime / numberOfIterations << " Sec\t " << DLP->actualMultiplicationTime / numberOfIterations << endl;

        cout << "\n Cheon Number Of Iterations :: " << numberOfIterations << " Done...\n";

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->teske3();
            teskeTime += DLP->getTimeByTeske();
        }
        teskeTime = teskeTime / numberOfIterations;

        if (whileLoopCnt == 0) {
            teske << "Number of Iterations :: " << numberOfIterations << endl;
            teske << "Number of Iterations of Walk :: " << DLP->getNumberOfIterations() << endl;
            teske << "\nr \tp^n \tTime Teske \t" << endl;
        }
        teske << r << "\t2^" << n << "\t" << teskeTime << " Sec" << endl;
        cout << "\n Teske Number Of Iterations :: " << numberOfIterations << " Done...\n";
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        DLP->~discreteLog();

        ++whileLoopCnt;
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}