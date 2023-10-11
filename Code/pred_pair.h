// pred_pair.h
#ifndef PRED_PAIR_H
#define PRED_PAIR_H
#include <stdint.h>
class PredPair
{
public:
    PredPair(): col_(0), value_(-1.0)
    {}
    PredPair(uint8_t col, double value): col_(col), value_(value)
    {}
    uint8_t col_;
    double value_;
    ~PredPair()
    {}
};
#endif