#include "tableCellGF2E.hpp"

tableCellGF2E::tableCellGF2E() {
    this->t = 0;
    this->p = 1;
}

void tableCellGF2E::printCellData() const {
    cout << "\n__________________________\n";
    cout << " Multiplier information ::";
    for (long int i = 0; i< this->numberOfElementsInMultiplierInformation; ++i) {
        cout << multiplierInformation[i] << " ";
    }
    cout << "\n groupElement  :: " << this->groupElement << endl;
    cout << " Summation ALPHA :: " << this->summationAlpha << endl;
    cout << " Summation BETA  :: " << this->summationBeta << endl;
    cout << " tag ::";
    for (int i = 0; i< this->n; ++i)
        cout << tag[i] << "\t";
    cout << "\n__________________________\n";
}

/**
 * @return Return's the tag for this cell in table Ml
 */
GF2E* tableCellGF2E::getTagFor() const {
    return this->tag;
}

void tableCellGF2E::setValues(long t, ZZ p, long long int l, ZZ n) {
    this->t = t;
    this->p = p;
    this->numberOfElementsInMultiplierInformation = l;
    this->n = n;

    multiplierInformation = new int[this->numberOfElementsInMultiplierInformation];

    tag = new GF2E[conv<long>(this->n)];
    //    for (int i = 0; i < this->n; ++i)
    //        tag[i].SetMaxLength(this->t);
}