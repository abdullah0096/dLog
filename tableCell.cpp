#include "tableCell.hpp"

void tableCell::printCellData() {
    cout << " Multiplier information ::";
    for (long int i = 0; i< this->l; ++i) {
        cout << multiplierInformation[i] << " ";
    }
    cout << "\n";
}

tableCell::tableCell() {
    this->t = 0;
    this->p = 1;
    //std::cout << "\n in here....\n";
}

void tableCell::setValues(long t, ZZ p, long long int l) {
    this->t = t;
    this->p = p;
    this->l = l;

    multiplierInformation = new int[this->l];
    tag = new int[this->t];
    ZZ_p::init(this->p);
}