#include "multiplierGF2E.hpp"

/**
 * Constructor accepts r as I/P i.e the size of set of multipliers also allocates
 * Memory for alpha, beta, i, and groupElements.
 * @param r : size of the set of Multipliers 
 * @TODO : A check for failure in memory allocation in case of large r.
 */
multiplierGF2E::multiplierGF2E(long r, const ZZ& p) {

    this->p = p;
    ZZ_p::init(this->p);

    this->r = r;
    alpha = new ZZ[this->r];
    beta = new ZZ[this->r];
    groupElement = new GF2E[this->r];
}

void multiplierGF2E::printMultiplier() const {

    for (long j = 0; j < r; ++j) {
        std::cout << "M[" << j << "] = " << this->groupElement[j] << "\t aplha :: " << this->alpha[j] << "\t beta :: " << this->beta[j] << std::endl;
    }
}