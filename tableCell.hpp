/* 
 * File:   tableCell.hpp
 * Author: Abdullah
 *
 * Created on October 21, 2014, 12:15 PM
 */

#ifndef TABLECELL_HPP
#define	TABLECELL_HPP
#include <iostream>

using namespace std;

class tableCell {
public:
    int *multiplierInformation; // Vector to store multiplier information
    //    NTL groupElement;
    long int summationAlpha; // summation of all alph's 
    long int summationBeta; // summation of all beta's 
    //    NTL tag;

public:

    tableCell() {

        std::cout << "\n in here.,,,\n";
    }

    void printCellData();


};

#endif	/* TABLECELL_HPP */

