#include "multiplier.hpp"

multiplier::multiplier() {
    this->alpha = 0;
    this->beta = 0;
    this->i = 0;
}

void multiplier::printMultiplier() {
    std::cout << "m[" << i << "]\talpha :: " << alpha << "\tbeta :: " << beta << std::endl;
}