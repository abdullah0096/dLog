#include "multiplier.hpp"

/**
 * Constructor accepts r as I/P i.e the size of set of multipliers also allocates
 * Memory for alpha, beta, i, and groupElements.
 * @param r : size of the set of Multipliers 
 * @TODO : A check for failure in memory allocation in case of large r.
 */
multiplier::multiplier(long r, const ZZ& p) {

    this->p = p;
    //    std::cout << "\n r :: " << r << "\t p :: " << p << std::endl;
    ZZ_p::init(this->p);

    this->r = r;
    alpha = new ZZ[this->r];
    beta = new ZZ[this->r];
    i = new long [this->r];
    groupElement = new ZZ_pX[this->r];
}

void multiplier::printMultiplier() {

    for (long j = 0; j < r; ++j) {
        std::cout << "M[" << this->i[j] << "] = " << this->groupElement[j] << "\t aplha :: " << this->alpha[j] << "\t beta :: " << this->beta[j] << std::endl;
    }
}

void multiplier::setAlpha(long i, long value) {
    this->alpha[i] = value;
}