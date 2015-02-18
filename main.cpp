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
#include <sstream>

#include "utility.hpp"

#include "discreteLog.hpp"
#include "discreteLogGF2.hpp"
#include "discreteLogGF2E.hpp"

using namespace NTL;
using namespace std;

long numberOfIterations = 100;
long numberOfIput = 2;

void ZZ_pX_W_MAT_W_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_ZZ_pX_W_MAT_W_S.txt");
    ofstream teske("exp2/teske_time_test2_ZZ_pX_W_MAT_W_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        ZZ_pX g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);
        fin>>irrdPoly;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
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
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void ZZ_pX_W_MAT_WO_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_ZZ_pX_W_MAT_WO_S.txt");
    ofstream teske("exp2/teske_time_test2_ZZ_pX_W_MAT_WO_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        ZZ_pX g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);
        fin >> irrdPoly;
        BuildIrred(irrdPoly, conv<long>(n));
        //        ZZ_pX id2;
        //        stringstream out;
        //        out << id2;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
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
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void ZZ_pX_WO_MAT_W_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_ZZ_pX_WO_MAT_W_S.txt");
    ofstream teske("exp2/teske_time_test2_ZZ_pX_WO_MAT_W_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        ZZ_pX g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);
        fin>>irrdPoly;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
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
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void ZZ_pX_WO_MAT_WO_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_ZZ_pX_WO_MAT_WO_S.txt");
    ofstream teske("exp2/teske_time_test2_ZZ_pX_WO_MAT_WO_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        ZZ_pX g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);
        ZZ_pX id2;
        fin >> irrdPoly;
        BuildIrred(irrdPoly, conv<long>(n));
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
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
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void GF2X_W_MAT_W_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_GF2X_W_MAT_W_S.txt");
    ofstream teske("exp2/teske_time_test2_GF2X_W_MAT_W_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        GF2X g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);
        fin>>irrdPoly;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLogGF2 *DLP;
        DLP = new discreteLogGF2(p, n, r, l, g, h, irrdPoly, t, orderOfG);

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->cheonDL2();
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
            DLP->teske2();
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
        DLP->~discreteLogGF2();

        ++whileLoopCnt;
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void GF2X_W_MAT_WO_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_GF2X_W_MAT_WO_S.txt");
    ofstream teske("exp2/teske_time_test2_GF2X_W_MAT_WO_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        GF2X g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);

        ZZ_pX id2;
        BuildIrred(id2, conv<long>(n));
        stringstream out;
        out << id2;
        fin >>irrdPoly;
        irrdPoly.zero();
        out>>irrdPoly;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLogGF2 *DLP;
        DLP = new discreteLogGF2(p, n, r, l, g, h, irrdPoly, t, orderOfG);

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->cheonDL2();
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
            DLP->teske2();
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
        DLP->~discreteLogGF2();

        ++whileLoopCnt;
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void GF2X_WO_MAT_W_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_GF2X_WO_MAT_W_S.txt");
    ofstream teske("exp2/teske_time_test2_GF2X_WO_MAT_W_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        GF2X g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);
        fin>>irrdPoly;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLogGF2 *DLP;
        DLP = new discreteLogGF2(p, n, r, l, g, h, irrdPoly, t, orderOfG);

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->cheonDL2();
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
            DLP->teske2();
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
        DLP->~discreteLogGF2();

        ++whileLoopCnt;
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void GF2X_WO_MAT_WO_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_GF2X_WO_MAT_WO_S.txt");
    ofstream teske("exp2/teske_time_test2_GF2X_WO_MAT_WO_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        GF2X g, h, irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);
        ZZ_pX id2;
        BuildIrred(id2, conv<long>(n));
        stringstream out;
        out << id2;
        fin >>irrdPoly;
        irrdPoly.zero();

        out>>irrdPoly;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLogGF2 *DLP;
        DLP = new discreteLogGF2(p, n, r, l, g, h, irrdPoly, t, orderOfG);

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->cheonDL2();
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
            DLP->teske2();
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
        DLP->~discreteLogGF2();

        ++whileLoopCnt;
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void GF2E_W_MAT_W_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_GF2E_W_MAT_W_S.txt");
    ofstream teske("exp2/teske_time_test2_GF2E_W_MAT_W_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;

        GF2X irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        ZZ_p::init(p);
        fin>>irrdPoly;
        GF2E::init(irrdPoly);
        GF2E g, h;
        fin >> g >> h;


        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLogGF2E *DLP;
        DLP = new discreteLogGF2E(p, n, r, l, g, h, irrdPoly, t, orderOfG);

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->cheonDL_Mat();
            cheonTime += DLP->getTimeByCheon();
        }

        cheonTime = cheonTime / numberOfIterations;
        if (whileLoopCnt == 0) {
            cheon << "Number of Iterations :: " << numberOfIterations << endl;
            cheon << "Number of Iterations of Walk :: " << DLP->getNumberOfIterations() << endl;
            cheon << "\nr l  t p^n \tTime Cheon \tTable Generation Time \t Gamma Time \tInner Prod Time \tTable Look-Up Time \tMiscellaneous Time \t Actual Multiplication" << endl;
        }

        cheon << std::fixed;
        cheon << r << " " << l << " " << trunc(log2(r)) << " 2^" << n << "\t" << cheonTime << " Sec\t   " << DLP->getTableGenerationTime() << " Sec\t\t "
                << DLP->cheon_gammaTime / numberOfIterations << " Sec\t " << DLP->innerProductTime / numberOfIterations
                << " Sec\t\t " << DLP->tableLookUpTime / numberOfIterations << " Sec\t\t"
                << DLP->cheon_miscellaneousTime / numberOfIterations << " Sec\t " << DLP->cheon_actualMultiplicationTime / numberOfIterations << endl;

        cout << "\n Cheon Number Of Iterations :: " << numberOfIterations << " Done...\n";

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->teske2();
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
        DLP->~discreteLogGF2E();

        ++whileLoopCnt;
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void GF2E_W_MAT_WO_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_GF2E_W_MAT_WO_S.txt");
    ofstream teske("exp2/teske_time_test2_GF2E_W_MAT_WO_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        GF2X irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);

        ZZ_pX id2;
        BuildIrred(id2, conv<long>(n));
        stringstream out;
        out << id2;
        fin >>irrdPoly;
        irrdPoly.zero();
        out>>irrdPoly;
        GF2E::init(irrdPoly);
        GF2E g, h;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLogGF2E *DLP;
        DLP = new discreteLogGF2E(p, n, r, l, g, h, irrdPoly, t, orderOfG);

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->cheonDL_Mat();
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
                << DLP->cheon_gammaTime / numberOfIterations << " Sec\t " << DLP->innerProductTime / numberOfIterations
                << " Sec\t\t " << DLP->tableLookUpTime / numberOfIterations << " Sec\t\t"
                << DLP->cheon_miscellaneousTime / numberOfIterations << " Sec\t " << DLP->cheon_actualMultiplicationTime / numberOfIterations << endl;

        cout << "\n Cheon Number Of Iterations :: " << numberOfIterations << " Done...\n";

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->teske2();
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
        DLP->~discreteLogGF2E();

        ++whileLoopCnt;
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void GF2E_WO_MAT_W_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_GF2E_WO_MAT_W_S.txt");
    ofstream teske("exp2/teske_time_test2_GF2E_WO_MAT_W_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        GF2X irrdPoly;
        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);

        fin>>irrdPoly;
        GF2E::init(irrdPoly);

        GF2E g, h;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLogGF2E *DLP;
        DLP = new discreteLogGF2E(p, n, r, l, g, h, irrdPoly, t, orderOfG);

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
                << DLP->cheon_gammaTime / numberOfIterations << " Sec\t " << DLP->innerProductTime / numberOfIterations
                << " Sec\t\t " << DLP->tableLookUpTime / numberOfIterations << " Sec\t\t"
                << DLP->cheon_miscellaneousTime / numberOfIterations << " Sec\t " << DLP->cheon_actualMultiplicationTime / numberOfIterations << endl;

        cout << "\n Cheon Number Of Iterations :: " << numberOfIterations << " Done...\n";

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->teske2();
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
        DLP->~discreteLogGF2E();

        ++whileLoopCnt;
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

void GF2E_WO_MAT_WO_S() {
    long double time = 0.0;
    long whileLoopCnt(0);
    ifstream fin("in2.txt");
    if (!fin) {
        cout << "\n ERROR in Main reading File in.txt...\n";
        exit(1);
    }
    ofstream cheon("exp2/cheon_time_test2_GF2E_WO_MAT_WO_S.txt");
    ofstream teske("exp2/teske_time_test2_GF2E_WO_MAT_WO_S.txt");

    while (!fin.eof()) {
        long r, l, t;
        ZZ p, n, orderOfG;
        GF2X irrdPoly;

        fin >> p >> n >> r >> orderOfG >> l >>t;

        ZZ_p::init(p);

        ZZ_pX id2;
        BuildIrred(id2, conv<long>(n));
        stringstream out;
        out << id2;
        fin >>irrdPoly;
        irrdPoly.zero();
        out>>irrdPoly;
        GF2E::init(irrdPoly);
        GF2E g, h;
        fin >> g >> h;

        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
                << "\t l :: " << l << "\t t ::" << t << "\t";
        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;

        long double cheonTime = 0;
        long double teskeTime = 0;
        discreteLogGF2E *DLP;
        DLP = new discreteLogGF2E(p, n, r, l, g, h, irrdPoly, t, orderOfG);

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
                << DLP->cheon_gammaTime / numberOfIterations << " Sec\t " << DLP->innerProductTime / numberOfIterations
                << " Sec\t\t " << DLP->tableLookUpTime / numberOfIterations << " Sec\t\t"
                << DLP->cheon_miscellaneousTime / numberOfIterations << " Sec\t " << DLP->cheon_actualMultiplicationTime / numberOfIterations << endl;

        cout << "\n Cheon Number Of Iterations :: " << numberOfIterations << " Done...\n";

        for (int i = 0; i < numberOfIterations; ++i) {
            DLP->teske2();
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
        DLP->~discreteLogGF2E();

        ++whileLoopCnt;
        if (whileLoopCnt >= numberOfIput) {
            fin.close();
            cheon.close();
            teske.close();
            break;
        }
    }//end::WHILE LOOP
    cheon.close();
    teske.close();
    fin.close();
}

int main(int argc, char** argv) {

    ZZ_pX_W_MAT_W_S();
    cout << "\n\n\n\t\t ZZ_pX_W_MAT_W_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    ZZ_pX_W_MAT_WO_S();
    cout << "\n\n\n\t\t ZZ_pX_W_MAT_WO_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    ZZ_pX_WO_MAT_W_S();
    cout << "\n\n\n\t\t ZZ_pX_WO_MAT_W_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    ZZ_pX_WO_MAT_WO_S();
    cout << "\n\n\n\t\t ZZ_pX_WO_MAT_WO_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    ;
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    ;
    GF2X_W_MAT_W_S();
    cout << "\n\n\n\t\t GF2X_W_MAT_W_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    GF2X_W_MAT_WO_S();
    cout << "\n\n\n\t\t GF2X_W_MAT_WO_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    GF2X_WO_MAT_W_S();
    cout << "\n\n\n\t\t GF2X_WO_MAT_W_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    GF2X_WO_MAT_WO_S();
    cout << "\n\n\n\t\t GF2X_WO_MAT_WO_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    ;
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    ;
    GF2E_W_MAT_W_S();
    cout << "\n\n\n\t\t GF2E_W_MAT_W_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    GF2E_W_MAT_WO_S();
    cout << "\n\n\n\t\t GF2E_W_MAT_WO_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    GF2E_WO_MAT_W_S();
    cout << "\n\n\n\t\t GF2E_WO_MAT_W_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    GF2E_WO_MAT_WO_S();
    cout << "\n\n\n\t\t GF2E_WO_MAT_WO_S() ... [DONE]";
    cout << "\n------------------------------------------------------------------\n\n\n";
    ;
    // <editor-fold defaultstate="collapsed" desc="GF2 Code ">
    //    long r, l, t;
    //    ZZ p, n, orderOfG;
    //    GF2X irrdPoly;
    //
    //    long numberOfIterations = 1;
    //
    //    ifstream fin("in.txt");
    //    if (!fin) {
    //        cout << "\n ERROR in Main reading File in.txt...\n";
    //        exit(1);
    //    }
    //    cout << "\n Reading Input from in.txt\n";
    //
    //    long r, l, t;
    //    ZZ p, n, orderOfG;
    //    GF2X irrdPoly;
    //    fin >> p >> n >> r >> orderOfG >> l >>t;
    //    ZZ_p::init(p);
    //
    //    ZZ_pX id2;
    //    BuildIrred(id2, conv<long>(n));
    //    stringstream out;
    //    out << id2;
    //    fin >>irrdPoly;
    //    out>>irrdPoly;
    //
    //    GF2E::init(irrdPoly);
    //    GF2E g, h;
    //
    //    fin >> g >> h;
    //
    //    long double time = 0;
    //
    //    discreteLogGF2E DLP(p, n, r, l, g, h, irrdPoly, t, orderOfG);
    //    DLP.printParameters();
    //    DLP.cheonDL3();
    //    cout << "\n Time :: " << DLP.getTimeByCheon() << endl;
    //    cout << "\n Cheon Time :: " << DLP.getTimeByCheon() << " Seconds." << endl;
    //    cout << "\nr l  t p^n \tTime Cheon \tTable Generation Time \t Gamma Time \tInner Prod Time \tTable Look-Up Time \tMiscellaneous Time \t Actual Multiplication" << endl;
    //    cout << r << " " << l << " " << trunc(log2(r)) << " 2^" << n << "\t" << DLP.getTimeByCheon()
    //            << " Sec\t   " << DLP.getTableGenerationTime() << " Sec\t\t "
    //            << DLP.cheon_gammaTime / numberOfIterations << " Sec\t " << DLP.innerProductTime / numberOfIterations
    //            << " Sec\t\t " << DLP.tableLookUpTime / numberOfIterations << " Sec\t\t"
    //            << DLP.cheon_miscellaneousTime / numberOfIterations << " Sec\t\t " << DLP.cheon_actualMultiplicationTime / numberOfIterations << endl;
    //    cout << r << " " << l << " " << trunc(log2(r)) << " 2^" << n << "\t" << DLP.getTimeByCheon()
    //            << " Sec\t   " << DLP.getTableGenerationTime() << " Sec\t\t "
    //            << ((DLP.cheon_gammaTime / numberOfIterations) / DLP.getTimeByCheon()) * 100 << " %  \t " << ((DLP.innerProductTime / numberOfIterations) / DLP.getTimeByCheon()) * 100
    //            << " %  \t\t " << ((DLP.tableLookUpTime / numberOfIterations) / DLP.getTimeByCheon()) * 100 << " %   \t\t"
    //            << ((DLP.cheon_miscellaneousTime / numberOfIterations) / DLP.getTimeByCheon()) * 100 << " %   \t\t " << ((DLP.cheon_actualMultiplicationTime / numberOfIterations) / DLP.getTimeByCheon()) * 100 << " %" << endl;
    //
    //    cout << "\n Collision Time :: " << DLP.collisionTime << endl;
    //    cout << "\n Walk Cnt Time :: " << DLP.cheon_walkCntTime << endl;
    //    cout << "\n-------------------------------------------------\n";
    //
    //    DLP.teske2();
    //
    //    cout << "\n Teske Time :: " << DLP.getTimeByTeske() << " Seconds." << endl;
    //    cout << "\n Gamma Time :: " << DLP.teske_gammaTime << "\t ( " << (DLP.teske_gammaTime / DLP.getTimeByTeske()) * 100 << " % )" << endl;
    //    cout << "\n Actual Mul :: " << DLP.teske_actualMultiplicationTime << "\t ( " << (DLP.teske_actualMultiplicationTime / DLP.getTimeByTeske())*100 << " % )" << endl;

    // </editor-fold>
    return 0;
}