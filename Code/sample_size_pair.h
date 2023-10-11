// sample_size_pair.h
#ifndef SAMPLE_SIZE_PAIR_H
#define SAMPLE_SIZE_PAIR_H
class SampleSizePair
{
public:
    SampleSizePair(): previous_(0), current_(0)
    {}
    SampleSizePair(int previous, int current): previous_(previous), current_(current)
    {}
    int previous_, current_;
    ~SampleSizePair()
    {}
};
#endif