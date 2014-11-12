#include "tableCell.hpp"

void tableCell::printCellData() {
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

tableCell::tableCell() {
    this->t = 0;
    this->p = 1;
    //std::cout << "\n in here....\n";
}

void tableCell::setValues(long t, ZZ p, long long int l, long long int n) {
    this->t = t;
    this->p = p;
    this->numberOfElementsInMultiplierInformation = l;
    this->n = n;

    multiplierInformation = new int[this->numberOfElementsInMultiplierInformation];
    ZZ_p::init(this->p);

    tag = new ZZ_pX[this->n];
    for (int i = 0; i < this->n; ++i)
        tag[i].SetMaxLength(this->t);
}

ZZ_pX* tableCell::getTagFor(){
    return this->tag;
}