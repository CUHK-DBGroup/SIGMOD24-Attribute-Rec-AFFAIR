// feature.h
#ifndef FEATURE_H
#define FEATURE_H
#include <stdlib.h> 
// const double small_double = 0.000000000001;
class Feature
{
public:
    Feature(): idx_(0), sample_size_(0), low_(0.0), up_(0.0), est_(0.0), avg_(0.0), alt_(0.0), diff_(1.0), base_low_(0.0), base_up_(0.0), to_update_(false)// , is_bound_(true)
    {}
    Feature(const uint32_t idx, const double init_value): idx_(idx), sample_size_(0), low_(init_value), up_(init_value), est_(init_value), avg_(init_value), alt_(init_value), diff_(1.0), to_update_(false)
    {}
    Feature(const uint32_t idx, const double low,  const double up, const double est, const uint32_t sample_size): idx_(idx), sample_size_(sample_size), low_(low), up_(up), est_(est), avg_((low_ + up_) / 2.0), alt_(up), diff_(up_ - low_), to_update_(false)
    {
        // avg_ = (low_ + up_) / 2.0;
        // diff_ = up_ - low_;
        // est_ = sample_est;
        // alt_ = up_;
    }
    uint32_t idx_, sample_size_;
    double low_, up_, est_, avg_, alt_, diff_, base_low_, base_up_, new_low_, new_up_;
    bool to_update_;
    void Update(const double low, const double up, const double est)
    {
        low_ = low;
        up_ = up;
        avg_ = (low_ + up_) / 2.0;
        diff_ = up_ - low_;
        est_ = est;
    }
    void UpdateWithSampleSize(const double low, const double up, const double diff, const double est, const uint32_t sample_size)
    {
        low_ = low;
        up_ = up;
        avg_ = (low_ + up_) / 2.0;
        diff_ = diff;
        est_ = est;
        alt_ = up_;
        sample_size_ = sample_size;
        to_update_ = false;
        // if (up_ > 0.0)
        // {   
        //     ratio_ = low_ / up_;
        // }
        // else
        // {
        //     ratio_ = 0.0;
        // }
    }
    void UpdateWithBaseline(const double low, const double up, const double sample_est, const double base_low, const double base_up, const double new_low, const double new_up)
    {
        low_ = low;
        up_ = up;
        // avg_ = (low_ + up_) / 2.0;
        diff_ = abs(up_ - low_);
        est_ = sample_est;
        base_low_ = base_low;
        base_up_ = base_up;
        new_low_ = new_low;
        new_up_ = new_up;
        // if (up_ > 0.0)
        // {   
        //     ratio_ = low_ / up_;
        // }
        // else
        // {
        //     ratio_ = 0.0;
        // }
    }
};
#endif