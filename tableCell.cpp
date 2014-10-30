#include "tableCell.hpp"

void tableCell::printCellData() {
    std::cout << "data";
}

tableCell::tableCell() {
    this->t = 0;
    this->p = 1;
    std::cout << "\n in here....\n";
}

void tableCell::setValues(long t, ZZ p) {
    this->t = t;
    this->p = p;
    ZZ_p::init(this->p);
}
