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

void InnerProduct1(ZZ_p& x, const Vec<ZZ_p>& a, const Vec<ZZ_p>& b) {
    long n = min(a.MaxLength(), b.MaxLength());
    //    cout << "\n n :: " << n << endl;
    long i;
    ZZ accum, t;

    accum = 0;
    for (i = 0; i < n; i++) {
        mul(t, rep(a[i]), rep(b[i]));
        add(accum, accum, t);
    }
    conv(x, accum);
    //    cout << "\n x :: " << x << endl;
}

void foo2() {
    ZZ_p::init(conv<ZZ>("564879132564879564132"));
    ZZ_p ans;
    vec_ZZ_p a, b;

    a.SetMaxLength(1024);
    b.SetMaxLength(1024);

    //    cout << "\n a :: " << a << endl;
    //    a.put(1, conv<ZZ_p>(1));

    SetSeed(conv<ZZ>("9343"));
    for (int i = 0; i < 1024; ++i) {
        a[i] = random_ZZ_p();
    }

    SetSeed(conv<ZZ>("3113"));
    for (int i = 0; i < 1024; ++i) {
        b[i] = random_ZZ_p();
    }

    //    for (int i = 0; i < 10; ++i) {
    //        cout << "a[" << i << "] :: " << a[i] << "\t b[" << i << "] :: " << b[i] << endl;
    //    }
    //    cout << "\n a.len() :: " << a.MaxLength() << "\t b.len() :: " << b.MaxLength() << endl;
    long n = min(a.MaxLength(), b.MaxLength());
    //    cout << "\n n :: " << n << endl;
    //    InnerProduct(ans, a, b);
    cout << "\n ans :: " << ans << endl;
}

int main(int argc, char** argv) {

    //    for (int i = 0; i < 1000; i++) {
    //        foo2();
    //    }

    //    InnerProduct1()
    foo();
    ;
    // <editor-fold defaultstate="collapsed" desc=" made fun foo ">
    //    long double time = 0.0;
    //    int end = 1;
    //
    //    for (long int i = 0; i < end; i++) {
    //        long r, l, t;
    //        ZZ p, n, orderOfG;
    //        ZZ_pX g, h, irrdPoly;
    //
    //        ifstream fin("in.txt");
    //        if (!fin) {
    //            cout << "\n ERROR in Main reading File in.txt...\n";
    //            exit(1);
    //        }
    //
    //        cout << "\n Reading Input from in.txt\n";
    //        fin >> p >> n >> r >> orderOfG >> l >>t;
    //
    //        ZZ_p::init(p);
    //
    //        fin >> g >> h>>irrdPoly;
    //        //        cout << "\n p :: " << p << "\t n :: " << n << "\t r :: " << r << "\t orderOfG :: " << orderOfG
    //        //                << "\t l :: " << l << "\t t ::" << t << "\t";
    //        //        cout << "\n g :: " << g << "\t h :: " << h << "\t irrdPoly :: " << irrdPoly << endl;
    //
    //        discreteLog DLP(p, n, r, l, g, h, irrdPoly, t, orderOfG);
    //        DLP.printParameters();
    //
    //        DLP.cheonDL3();
    //        cout << "\n Time By Cheon2 ::" << DLP.getTimeByCheon() << " Seconds..." << endl;
    //
    //        DLP.teske3();
    //        cout << "\n Time By teske ::" << DLP.getTimeByTeske() << " Seconds..." << endl;
    //
    //        fin.close();
    //    }
    // </editor-fold>
    ;
    // <editor-fold defaultstate="collapsed" desc="GF2 Code ">
    //    long r, l, t;
    //    ZZ p, n, orderOfG;
    //    GF2X g, h, irrdPoly;
    //
    //    ifstream fin("in.txt");
    //    if (!fin) {
    //        cout << "\n ERROR in Main reading File in.txt...\n";
    //        exit(1);
    //    }
    //
    //    cout << "\n Reading Input from in.txt\n";
    //    fin >> p >> n >> r >> orderOfG >> l >>t;
    //
    //    ZZ_p::init(p);
    //    fin >> g >> h>>irrdPoly;
    //
    //    discreteLogGF2 DLP(p, n, r, l, g, h, irrdPoly, t, orderOfG);
    //    DLP.printParameters();
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
    ofstream cheon("cheon_time_test.txt");
    ofstream teske("teske_time_test.txt");

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