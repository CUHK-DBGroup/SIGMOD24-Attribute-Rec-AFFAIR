// config.h
#ifndef	CONFIG_H
#define CONFIG_H
#include <string>
#include "pred_pair.h"
class Config 
{
public:
    std::string action_ = "", prefix_ = "/home/Dataset/", dataset_ = "cdc82_01", metric_ = "ks";
    double empirical_threshold_ = 0.0, abs_error_ = 0.0, fail_pr_;
    uint32_t topk_ = 4, initial_size_ = 1024, query_num_ = 1, batch_size_ = 10000000;
    std::pair<double, double> predicate1_range_ratio_{0.0, 0.5}, predicate2_range_ratio_{0.5, 1.0};
    std::pair<uint32_t, uint32_t> predicate1_range_row_{0, 0}, predicate2_range_row_{0, 0}, predicate_size_{0, 0};
    Config()
    {}
};
#endif