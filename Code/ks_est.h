// ks_est.h
#ifndef KS_EST_H
#define METRIC_EST_H
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>
#include <string>
#include "feature.h"
#include "sample_size_pair.h"

#define likely(x) __builtin_expect(!!(x), 1) 
#define unlikely(x) __builtin_expect(!!(x), 0)
class MetricEst
{
private:
    double ExpectedDiffCd(const double& integral_tail_pr, std::vector<uint32_t>& sample_hist, const uint32_t sample_size, const uint32_t support_size)
    {
        const double init_fair_pr = 0.5;
        double expected_diff = 0.0, single_fail_pr = init_fair_pr;
        while (single_fail_pr >= integral_tail_pr)
        {
            // double max_diff = 1.0;
            // double a = log(2.0 / single_fail_pr);
            double max_diff = 1.0, a = log(support_size * 2.0 / single_fail_pr);
            for (uint32_t i = 0; i < sample_hist.size(); i++)
            {
                double est_x = 1.0 * sample_hist[i] / sample_size, 
                lower_bound = std::max(0.0, (est_x + 2.0 * a / 3.0 / sample_size - sqrt(4.0 * a * a / 3.0 / sample_size / sample_size * (1.0 / 3.0 + est_x) + 2.0 * a / sample_size * est_x * (1.0 - est_x))) / (2.0 * a / sample_size + 1.0)),
                upper_bound = std::min(1.0, est_x + 1.0 * a / sample_size + sqrt(2.0 * a / sample_size * (est_x + a / 2.0 / sample_size))),
                // curr_diff = upper_bound - lower_bound;
                curr_diff = std::max(est_x - lower_bound, upper_bound - est_x);
                // printf("new\n");
                if (curr_diff < max_diff)
                {
                    max_diff = curr_diff;
                }
            }
            // printf("max_diff is %f\n", max_diff);
            expected_diff += max_diff * single_fail_pr;
            single_fail_pr /= 2.0;
        }
        expected_diff += 1.0 * single_fail_pr;
        return expected_diff;
    }

    double ExpectedDiffKs(const double& integral_tail_pr, std::vector<uint32_t>& sample_hist, const uint32_t sample_size, const uint32_t support_size)
    {
        const double init_fair_pr = 0.5;
        double expected_diff = 0.0, single_fail_pr = init_fair_pr;
        while (single_fail_pr >= integral_tail_pr)
        {
            // double max_diff = 1.0, a = log(2.0 / single_fail_pr), est_x = 0.0;
            double max_diff = 1.0, a = log(support_size * 2.0 / single_fail_pr), est_x = 0.0;
            // printf("new 1\n");
            for (uint32_t i = 0; i < sample_hist.size(); i++)
            {
                est_x += 1.0 * sample_hist[i] / sample_size;
                double lower_bound = std::max(0.0, (est_x + 2.0 * a / 3.0 / sample_size - sqrt(4.0 * a * a / 3.0 / sample_size / sample_size * (1.0 / 3.0 + est_x) + 2.0 * a / sample_size * est_x * (1.0 - est_x))) / (2.0 * a / sample_size + 1.0)),
                upper_bound = std::min(1.0, est_x + 1.0 * a / sample_size + sqrt(2.0 * a / sample_size * (est_x + a / 2.0 / sample_size))),
                // curr_diff = upper_bound - lower_bound;
                curr_diff = std::max(est_x - lower_bound, upper_bound - est_x);
                if (curr_diff < max_diff)
                {
                    max_diff = curr_diff;
                }
            }
            // printf("max_diff is %f\n", max_diff);
            expected_diff += max_diff * single_fail_pr;
            single_fail_pr /= 2.0;
        }
        expected_diff += 1.0 * single_fail_pr;
        return expected_diff;
    }

    double ExpectedDiffEMD(const double& integral_tail_pr, std::vector<uint32_t>& sample_hist, const uint32_t sample_size, const uint32_t support_size)
    {
        const double init_fair_pr = 0.5;
        double expected_diff = 0.0, single_fail_pr = init_fair_pr;
        while (single_fail_pr >= integral_tail_pr)
        {
            // double max_diff = 1.0, a = log(2.0 / single_fail_pr), est_x = 0.0;
            double sum_diff = 0.0, a = log(support_size * 2.0 / single_fail_pr);
            // printf("new 1\n");
            for (uint32_t i = 0; i < sample_hist.size(); i++)
            {
                double est_x = 1.0 * sample_hist[i] / sample_size, 
                lower_bound = std::max(0.0, (est_x + 2.0 * a / 3.0 / sample_size - sqrt(4.0 * a * a / 3.0 / sample_size / sample_size * (1.0 / 3.0 + est_x) + 2.0 * a / sample_size * est_x * (1.0 - est_x))) / (2.0 * a / sample_size + 1.0)),
                upper_bound = std::min(1.0, est_x + 1.0 * a / sample_size + sqrt(2.0 * a / sample_size * (est_x + a / 2.0 / sample_size))),
                // curr_diff = upper_bound - lower_bound;
                curr_diff = std::max(est_x - lower_bound, upper_bound - est_x);
                sum_diff += curr_diff;
            }
            // printf("max_diff is %f\n", max_diff);
            expected_diff += sum_diff / 2.0 * single_fail_pr;
            single_fail_pr /= 2.0;
        }
        expected_diff += 1.0 * single_fail_pr;
        return expected_diff;
    }

    double ExpectedDiffEuclidean(const double& integral_tail_pr, std::vector<uint32_t>& sample_hist, const uint32_t sample_size, const uint32_t support_size)
    {
        const double init_fair_pr = 0.5;
        double expected_diff = 0.0, single_fail_pr = init_fair_pr;
        while (single_fail_pr >= integral_tail_pr)
        {
            // double max_diff = 1.0, a = log(2.0 / single_fail_pr), est_x = 0.0;
            double sum_squared_diff = 0.0, a = log(support_size * 2.0 / single_fail_pr);
            // printf("new 1\n");
            for (uint32_t i = 0; i < sample_hist.size(); i++)
            {
                double est_x = 1.0 * sample_hist[i] / sample_size, 
                lower_bound = std::max(0.0, (est_x + 2.0 * a / 3.0 / sample_size - sqrt(4.0 * a * a / 3.0 / sample_size / sample_size * (1.0 / 3.0 + est_x) + 2.0 * a / sample_size * est_x * (1.0 - est_x))) / (2.0 * a / sample_size + 1.0)),
                upper_bound = std::min(1.0, est_x + 1.0 * a / sample_size + sqrt(2.0 * a / sample_size * (est_x + a / 2.0 / sample_size))),
                // curr_diff = upper_bound - lower_bound;
                curr_diff = std::max(est_x - lower_bound, upper_bound - est_x);
                sum_squared_diff += curr_diff * curr_diff;
            }
            // printf("max_diff is %f\n", max_diff);
            expected_diff += sqrt(sum_squared_diff / 2.0) * single_fail_pr;
            single_fail_pr /= 2.0;
        }
        expected_diff += 1.0 * single_fail_pr;
        return expected_diff;
    }

    double ExpectedDiffKsNew(const double& integral_tail_pr, std::vector<uint32_t>& sample_hist, const uint32_t sample_size, const uint32_t total_size)
    {
        const double init_fair_pr = 0.5;
        double expected_diff = 0.0, single_fail_pr = init_fair_pr, rho, kappa;
        if (sample_size <= total_size / 2) {
            rho = 1.0 - 1.0 * (sample_size - 1.0) / total_size;
            kappa = 4.0 / 3.0 + sqrt((1.0 * sample_size / total_size) / (total_size / (sample_size - 1.0) - 1.0));    
        }
        else {
            rho = (1.0 - 1.0 * sample_size / total_size) * (1.0 + 1.0 / sample_size);
            kappa = 4.0 / 3.0 + sqrt((1.0 - 1.0 * sample_size / total_size) * (total_size / (sample_size + 1.0) - 1.0));
        }
        while (single_fail_pr >= integral_tail_pr)
        {
            double max_diff = 1.0, alpha = 2.0 * rho * log(4.0 / single_fail_pr) / sample_size, beta = kappa * log(4.0 / single_fail_pr) / sample_size, est_x = 0.0;
            for (uint32_t i = 0; i < sample_hist.size(); i++)
            {
                est_x += 1.0 * sample_hist[i] / sample_size;
                double lower_bias = std::max(0.0, (2.0 * beta - alpha + 2.0 * alpha * est_x + sqrt((2.0 * beta - 2.0 * est_x - alpha) * (2.0 * beta - 2.0 * est_x - alpha) - 4.0 * (1.0 + alpha) * (beta * beta - 2.0 * beta * est_x + est_x * est_x))) / (2.0 * (1.0 + alpha))),
                upper_bias = std::max(0.0, (2.0 * beta + alpha + sqrt((2.0 * beta + 2.0 * est_x + alpha) * (2.0 * beta + 2.0 * est_x + alpha) - 4.0 * (1.0 + alpha) * (beta * beta + 2.0 * beta * est_x + est_x * est_x)) - 2.0 * alpha * est_x) / (2.0 * (1.0 + alpha))),
                curr_diff = std::min(1.0, std::max(lower_bias, upper_bias));
                if (curr_diff < max_diff)
                {
                    max_diff = curr_diff;
                }
            }
            // printf("max_diff is %f\n", max_diff);
            expected_diff += max_diff * single_fail_pr;
            single_fail_pr /= 2.0;
        }
        expected_diff += 1.0 * single_fail_pr;
        return expected_diff;
    }

    void UpdateSampleHist(const std::vector<uint32_t>& in_vec, std::vector<uint32_t>& in_hist, const std::vector<uint32_t>& predicate_sample_index) 
    {
        for (uint32_t i = 0; i < predicate_sample_index.size(); i++) 
        {
            in_hist[in_vec[predicate_sample_index[i]]] ++;
        }  
    }

    void UpdateSampleHist2(const std::vector<uint32_t>& in_vec, std::vector<uint32_t>& in_hist, const std::vector<uint32_t>& predicate_sample_index, const std::pair<uint32_t, uint32_t>& predicate_sample_size_range) 
    {
        for (uint32_t i = predicate_sample_size_range.first; i < predicate_sample_size_range.second; i++) 
        {
            in_hist[in_vec[predicate_sample_index[i]]] ++;
        }  
    }

    double SingleRelDiffIntegral(const double expected_bound, const uint32_t in_vec_size, const uint32_t sample_size, double fail_pr)
    {
        return expected_bound + sqrt(1.0 * (in_vec_size - sample_size) / (1.0 * sample_size * in_vec_size)) * sqrt(2.0 * log(2.0 / fail_pr));
    }
    
    double SingleRelDiffIntegralNew(const double expected_bound, const uint32_t in_vec_size, const uint32_t sample_size, double fail_pr)
    {
        return expected_bound + sqrt(1.0 * (in_vec_size - sample_size) / 4.0 / (1.0 * sample_size * (in_vec_size - 0.5) * (1.0 - 1.0 / std::max(sample_size, in_vec_size - sample_size)))) * sqrt(2.0 * log(2.0 / fail_pr));
    }
    
	int num_element_;
	double fail_pr_;
    std::string mode_;
    Feature variable_, variable_square_;
public:
	MetricEst(): num_element_(0), fail_pr_(0.25), low_(0.0), up_(0.0), est_(0.0) {}
    MetricEst(const std::string& mode, const Config& config, const std::vector<uint32_t>& in_vec, const std::pair<std::vector<uint32_t>, std::vector<uint32_t>>& predicate_sample_index_pair, const std::pair<uint32_t, uint32_t>& predicate_sample_size, const uint32_t support_size, const double fail_pr, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>& sample_hist_pair, uint32_t idx): 
        num_element_(in_vec.size()), fail_pr_(fail_pr), low_(0.0), up_(0.0), est_(0.0) 
    {
        // double rel_diff_integral = 0.0, rel_diff_baseline = 0.0, rel_diff_new = 0.0;
        UpdateSampleHist(in_vec, sample_hist_pair.first, predicate_sample_index_pair.first);
        UpdateSampleHist(in_vec, sample_hist_pair.second, predicate_sample_index_pair.second);
        // printf("pair size %d %d predicate sample size %d %d", predicate_sample_index_pair.first.size(), predicate_sample_index_pair.second.size(), predicate_sample_size.first, predicate_sample_size.second);
        double rel_diff = 0.0;
        if (mode == "approx")
        {
            // double rel_diff_integral = 0.0, rel_diff_baseline = 0.0, rel_diff_new = 0.0;
            std::pair<double, double> expected_bound_pair(0.0, 0.0);
            // printf("hist1 & hist2\n");

            if (config.metric_ == "ks") 
            {
                // std::pair<uint32_t, uint32_t> cum_hist_pair(0, 0), total_size_pair(in_vec1.size(), in_vec2.size()); 
                std::pair<uint32_t, uint32_t> cum_hist_pair(0, 0);
                for (int i = 0; i < support_size; i++) 
                {
                    cum_hist_pair.first += sample_hist_pair.first[i];
                    cum_hist_pair.second += sample_hist_pair.second[i];
                    double curr_diff_pr = abs(1.0 * cum_hist_pair.first / predicate_sample_size.first - 1.0 * cum_hist_pair.second / predicate_sample_size.second);
                    if (est_ < curr_diff_pr) 
                    {
                        est_ = curr_diff_pr;
                    }         
                }
                // ExpectedBoundIntegralKs(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                expected_bound_pair.first = ExpectedDiffKs(1.0 / config.predicate_size_.first, sample_hist_pair.first, predicate_sample_size.first, support_size);
                expected_bound_pair.second = ExpectedDiffKs(1.0 / config.predicate_size_.second, sample_hist_pair.second, predicate_sample_size.second, support_size);
            }
            else if (config.metric_ == "cd") 
            {
                for (int i = 0; i < support_size; i++) 
                {
                    // printf("%d, %d\n", sample_hist_pair.first[i], sample_hist_pair.second[i]);
                    double curr_diff_pr = abs(1.0 * sample_hist_pair.first[i] / predicate_sample_size.first - 1.0 * sample_hist_pair.second[i] / predicate_sample_size.second);
                    if (est_ < curr_diff_pr) 
                    {
                        est_ = curr_diff_pr;
                    }         
                }
                // ExpectedBoundIntegralCd(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                expected_bound_pair.first = ExpectedDiffCd(1.0 / config.predicate_size_.first, sample_hist_pair.first, predicate_sample_size.first, support_size);
                expected_bound_pair.second = ExpectedDiffCd(1.0 / config.predicate_size_.second, sample_hist_pair.second, predicate_sample_size.second, support_size);
            }
            else if (config.metric_ == "emd") 
            {
                double tmp_est = 0.0;
                for (int i = 0; i < support_size; i++) 
                {
                    double curr_diff_pr = abs(1.0 * sample_hist_pair.first[i] / predicate_sample_size.first - 1.0 * sample_hist_pair.second[i] / predicate_sample_size.second);
                    tmp_est += curr_diff_pr;
                }
                est_ = tmp_est / 2.0;
                // ExpectedBoundIntegralKs(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                expected_bound_pair.first = ExpectedDiffEMD(1.0 / config.predicate_size_.first, sample_hist_pair.first, predicate_sample_size.first, support_size);
                expected_bound_pair.second = ExpectedDiffEMD(1.0 / config.predicate_size_.second, sample_hist_pair.second, predicate_sample_size.second, support_size);
            }
            else if (config.metric_ == "euclidean") 
            {
                double tmp_est = 0.0;
                for (int i = 0; i < support_size; i++) 
                {
                    double curr_diff_pr = abs(1.0 * sample_hist_pair.first[i] / predicate_sample_size.first - 1.0 * sample_hist_pair.second[i] / predicate_sample_size.second);
                    tmp_est += curr_diff_pr * curr_diff_pr;
                }
                est_ = sqrt(tmp_est / 2.0);
                // ExpectedBoundIntegralKs(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                expected_bound_pair.first = ExpectedDiffEuclidean(1.0 / config.predicate_size_.first, sample_hist_pair.first, predicate_sample_size.first, support_size);
                expected_bound_pair.second = ExpectedDiffEuclidean(1.0 / config.predicate_size_.second, sample_hist_pair.second, predicate_sample_size.second, support_size);
            }
            else
            {
                printf("metric not recognize\n");
            }
            // printf("current sample_size is %d est_ks_test is %f\n", sample_pair.second, est_);
            // printf("expected bound first %f second %f\n", expected_bound_pair.first, expected_bound_pair.second);
            // rel_diff_integral += SingleRelDiffIntegral(expected_bound_pair.first, in_vec1.size(), sample_pair.second, fail_pr_);
            // rel_diff_integral += SingleRelDiffIntegral(expected_bound_pair.second, in_vec2.size(), sample_pair.second, fail_pr_);
            // rel_diff_baseline += 1.0 / 2.0 * sqrt((support_size - 1.0) * ((in_vec1.size() - std::min(sample_pair.second, uint32_t(in_vec1.size()))) / (1.0 * std::min(sample_pair.second, uint32_t(in_vec1.size())) * (in_vec1.size() - 1.0)))) + sqrt(1.0 * (in_vec1.size() - std::min(sample_pair.second, uint32_t(in_vec1.size()))) / (1.0 * std::min(sample_pair.second, uint32_t(in_vec1.size())) * in_vec1.size())) * sqrt(2.0 * log(2.0 / fail_pr_));
            // rel_diff_baseline += 1.0 / 2.0 * sqrt((support_size - 1.0) * ((in_vec2.size() - std::min(sample_pair.second, uint32_t(in_vec2.size()))) / (1.0 * std::min(sample_pair.second, uint32_t(in_vec2.size())) * (in_vec2.size() - 1.0)))) + sqrt(1.0 * (in_vec2.size() - std::min(sample_pair.second, uint32_t(in_vec2.size()))) / (1.0 * std::min(sample_pair.second, uint32_t(in_vec2.size())) * in_vec2.size())) * sqrt(2.0 * log(2.0 / fail_pr_));
            rel_diff += SingleRelDiffIntegralNew(expected_bound_pair.first, config.predicate_size_.first, predicate_sample_size.first, fail_pr_);
            rel_diff += SingleRelDiffIntegralNew(expected_bound_pair.second, config.predicate_size_.second, predicate_sample_size.second, fail_pr_);
        }
        else if (mode == "baseline")
        {
            if (config.metric_ == "emd")
            {
                double tmp_est = 0.0, 
                h1 = 1.0 / 2.0 * sqrt((support_size - 1.0) * ((config.predicate_size_.first - predicate_sample_size.first) / (1.0 * predicate_sample_size.first * (config.predicate_size_.first - 1.0)))), 
                h2 = 1.0 / 2.0 * sqrt((support_size - 1.0) * ((config.predicate_size_.second - predicate_sample_size.second) / (1.0 * predicate_sample_size.second * (config.predicate_size_.second - 1.0)))),
                g1 = sqrt(1.0 * (config.predicate_size_.first - predicate_sample_size.first) / (1.0 * predicate_sample_size.first * config.predicate_size_.first)),
                g2 = sqrt(1.0 * (config.predicate_size_.second - predicate_sample_size.second) / (1.0 * predicate_sample_size.second * config.predicate_size_.second));
                for (int i = 0; i < support_size; i++) 
                {
                    double curr_diff_pr = abs(1.0 * sample_hist_pair.first[i] / predicate_sample_size.first - 1.0 * sample_hist_pair.second[i] / predicate_sample_size.second);
                    tmp_est += curr_diff_pr;
                }
                est_ = tmp_est / 2.0;
                // rel_diff += h1 + g1 * sqrt(2.0 * log(sqrt(g1 + g2) / sqrt(g1) / fail_pr_));
                // rel_diff += h2 + g2 * sqrt(2.0 * log(sqrt(g1 + g2) / sqrt(g2) / fail_pr_));

                if (g1 > 0) 
                {
                    rel_diff += h1 + g1 * sqrt(2.0 * log(sqrt(g1 + g2) / sqrt(g1) / fail_pr_));
                }
                else if (g2 > 0)
                {
                    rel_diff += h2 + g2 * sqrt(2.0 * log(sqrt(g1 + g2) / sqrt(g2) / fail_pr_));
                }
                // rel_diff += 1.0 / 2.0 * sqrt((support_size - 1.0) * ((in_vec1.size() - std::min(sample_pair.second, uint32_t(in_vec1.size()))) / (1.0 * std::min(sample_pair.second, uint32_t(in_vec1.size())) * (in_vec1.size() - 1.0)))) + sqrt(1.0 * (in_vec1.size() - std::min(sample_pair.second, uint32_t(in_vec1.size()))) / (1.0 * std::min(sample_pair.second, uint32_t(in_vec1.size())) * in_vec1.size())) * sqrt(2.0 * log(2.0 / fail_pr_));
                // rel_diff += 1.0 / 2.0 * sqrt((support_size - 1.0) * ((in_vec2.size() - std::min(sample_pair.second, uint32_t(in_vec2.size()))) / (1.0 * std::min(sample_pair.second, uint32_t(in_vec2.size())) * (in_vec2.size() - 1.0)))) + sqrt(1.0 * (in_vec2.size() - std::min(sample_pair.second, uint32_t(in_vec2.size()))) / (1.0 * std::min(sample_pair.second, uint32_t(in_vec2.size())) * in_vec2.size())) * sqrt(2.0 * log(2.0 / fail_pr_));
                // rel_diff += SingleRelDiffIntegralNew(expected_bound_pair.first, in_vec1.size(), sample_pair.second, fail_pr_);
            }
            else if (config.metric_ == "euclidean") 
            {
                double tmp_est = 0.0,
                h1 = sqrt((support_size - 1.0) * ((config.predicate_size_.first - predicate_sample_size.first) / (2.0 * support_size * predicate_sample_size.first * (config.predicate_size_.first - 1.0)))), 
                h2 = sqrt((support_size - 1.0) * ((config.predicate_size_.second - predicate_sample_size.second) / (2.0 * support_size * predicate_sample_size.second * (config.predicate_size_.second - 1.0)))),
                g1 = sqrt(1.0 * (config.predicate_size_.first - predicate_sample_size.first) / (1.0 * predicate_sample_size.first * config.predicate_size_.first)),
                g2 = sqrt(1.0 * (config.predicate_size_.second - predicate_sample_size.second) / (1.0 * predicate_sample_size.second * config.predicate_size_.second));
                for (int i = 0; i < support_size; i++) 
                {
                    double curr_diff_pr = abs(1.0 * sample_hist_pair.first[i] / predicate_sample_size.first - 1.0 * sample_hist_pair.second[i] / predicate_sample_size.second);
                    tmp_est += curr_diff_pr * curr_diff_pr;
                }
                est_ = sqrt(tmp_est / 2.0);
                // rel_diff += h1 + g1 * sqrt(2.0 * log(sqrt(g1 + g2) / sqrt(g1) / fail_pr_));
                // rel_diff += h2 + g2 * sqrt(2.0 * log(sqrt(g1 + g2) / sqrt(g2) / fail_pr_));

                if (g1 > 0) 
                {
                    rel_diff += h1 + g1 * sqrt(2.0 * log(sqrt(g1 + g2) / sqrt(g1) / fail_pr_));
                }
                else if (g2 > 0)
                {
                    rel_diff += h2 + g2 * sqrt(2.0 * log(sqrt(g1 + g2) / sqrt(g2) / fail_pr_));
                }
                // ExpectedBoundIntegralKs(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                // expected_bound_pair.first = ExpectedDiffEuclidean(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair.first, sample_pair.second, support_size);
                // expected_bound_pair.second = ExpectedDiffEuclidean(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair.second, sample_pair.second, support_size);
            }
            else
            {
                printf("metric not recognize\n");
            }
        }
        else
        {
            printf("unrecognized mode\n");
        }
        
        if (predicate_sample_size.first >= config.predicate_size_.first || predicate_sample_size.second >= config.predicate_size_.second)
        {
            rel_diff = 0.0;
        }
        if (idx == 0 || idx == 35) 
        {
            printf("idx %d low_ %f up_ %f est_ %f diff_ %f\n", idx, low_, up_, est_, diff_);
        }
        diff_ = rel_diff;
        low_ = std::max(0.0, est_ - rel_diff);
        up_ = std::min(1.0, est_ + rel_diff);
        est_ = (est_ + rel_diff) > 1.0? (low_ + up_) / 2.0: est_;
    }
    MetricEst(const std::string& mode, const Config& config, const std::vector<uint32_t>& in_vec, const std::pair<std::vector<uint32_t>, std::vector<uint32_t>>& predicate_sample_index_pair, const std::pair<std::pair<uint32_t, uint32_t>, std::pair<uint32_t, uint32_t>>& predicate_sample_size_range, const uint32_t support_size, const double fail_pr, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>& sample_hist_pair, uint32_t idx): 
        num_element_(in_vec.size()), fail_pr_(fail_pr), low_(0.0), up_(0.0), est_(0.0) 
    {
        // double rel_diff_integral = 0.0, rel_diff_baseline = 0.0, rel_diff_new = 0.0;
        UpdateSampleHist2(in_vec, sample_hist_pair.first, predicate_sample_index_pair.first, predicate_sample_size_range.first);
        UpdateSampleHist2(in_vec, sample_hist_pair.second, predicate_sample_index_pair.second, predicate_sample_size_range.second);
        // printf("pair size %d %d predicate sample size %d %d", predicate_sample_index_pair.first.size(), predicate_sample_index_pair.second.size(), predicate_sample_size.first, predicate_sample_size.second);
        double rel_diff = 0.0;
        if (mode == "approx")
        {
            // double rel_diff_integral = 0.0, rel_diff_baseline = 0.0, rel_diff_new = 0.0;
            std::pair<double, double> expected_bound_pair(0.0, 0.0);
            // printf("hist1 & hist2\n");

            if (config.metric_ == "ks") 
            {
                // std::pair<uint32_t, uint32_t> cum_hist_pair(0, 0), total_size_pair(in_vec1.size(), in_vec2.size()); 
                std::pair<uint32_t, uint32_t> cum_hist_pair(0, 0);
                for (int i = 0; i < support_size; i++) 
                {
                    cum_hist_pair.first += sample_hist_pair.first[i];
                    cum_hist_pair.second += sample_hist_pair.second[i];
                    double curr_diff_pr = abs(1.0 * cum_hist_pair.first / predicate_sample_size_range.first.second - 1.0 * cum_hist_pair.second / predicate_sample_size_range.second.second);
                    if (est_ < curr_diff_pr) 
                    {
                        est_ = curr_diff_pr;
                    }         
                }
                // ExpectedBoundIntegralKs(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                expected_bound_pair.first = ExpectedDiffKs(1.0 / config.predicate_size_.first, sample_hist_pair.first, predicate_sample_size_range.first.second, support_size);
                expected_bound_pair.second = ExpectedDiffKs(1.0 / config.predicate_size_.second, sample_hist_pair.second, predicate_sample_size_range.second.second, support_size);
            }
            else if (config.metric_ == "cd") 
            {
                for (int i = 0; i < support_size; i++) 
                {
                    // printf("%d, %d\n", sample_hist_pair.first[i], sample_hist_pair.second[i]);
                    double curr_diff_pr = abs(1.0 * sample_hist_pair.first[i] / predicate_sample_size_range.first.second - 1.0 * sample_hist_pair.second[i] / predicate_sample_size_range.second.second);
                    if (est_ < curr_diff_pr) 
                    {
                        est_ = curr_diff_pr;
                    }         
                }
                // ExpectedBoundIntegralCd(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                expected_bound_pair.first = ExpectedDiffCd(1.0 / config.predicate_size_.first, sample_hist_pair.first, predicate_sample_size_range.first.second, support_size);
                expected_bound_pair.second = ExpectedDiffCd(1.0 / config.predicate_size_.second, sample_hist_pair.second, predicate_sample_size_range.second.second, support_size);
            }
            else if (config.metric_ == "emd") 
            {
                double tmp_est = 0.0;
                for (int i = 0; i < support_size; i++) 
                {
                    double curr_diff_pr = abs(1.0 * sample_hist_pair.first[i] / predicate_sample_size_range.first.second - 1.0 * sample_hist_pair.second[i] / predicate_sample_size_range.second.second);
                    tmp_est += curr_diff_pr;
                }
                est_ = tmp_est / 2.0;
                // ExpectedBoundIntegralKs(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                expected_bound_pair.first = ExpectedDiffEMD(1.0 / config.predicate_size_.first, sample_hist_pair.first, predicate_sample_size_range.first.second, support_size);
                expected_bound_pair.second = ExpectedDiffEMD(1.0 / config.predicate_size_.second, sample_hist_pair.second, predicate_sample_size_range.second.second, support_size);
            }
            else if (config.metric_ == "euclidean") 
            {
                double tmp_est = 0.0;
                for (int i = 0; i < support_size; i++) 
                {
                    double curr_diff_pr = abs(1.0 * sample_hist_pair.first[i] / predicate_sample_size_range.first.second - 1.0 * sample_hist_pair.second[i] / predicate_sample_size_range.second.second);
                    tmp_est += curr_diff_pr * curr_diff_pr;
                }
                est_ = sqrt(tmp_est / 2.0);
                // ExpectedBoundIntegralKs(1.0 / std::max(in_vec1.size(), in_vec2.size()), sample_hist_pair, sample_pair.second, expected_bound_pair, support_size);
                expected_bound_pair.first = ExpectedDiffEuclidean(1.0 / config.predicate_size_.first, sample_hist_pair.first, predicate_sample_size_range.first.second, support_size);
                expected_bound_pair.second = ExpectedDiffEuclidean(1.0 / config.predicate_size_.second, sample_hist_pair.second, predicate_sample_size_range.second.second, support_size);
            }
            else
            {
                printf("metric not recognize\n");
            }
            
            rel_diff += SingleRelDiffIntegralNew(expected_bound_pair.first, config.predicate_size_.first, predicate_sample_size_range.first.second, fail_pr_);
            rel_diff += SingleRelDiffIntegralNew(expected_bound_pair.second, config.predicate_size_.second, predicate_sample_size_range.second.second, fail_pr_);
        }
        else
        {
            printf("unrecognized mode\n");
        }
        
        if (predicate_sample_size_range.first.second >= config.predicate_size_.first || predicate_sample_size_range.second.second >= config.predicate_size_.second)
        {
            rel_diff = 0.0;
        }
        diff_ = rel_diff;
        low_ = std::max(0.0, est_ - rel_diff);
        up_ = std::min(1.0, est_ + rel_diff);
        est_ = (est_ + rel_diff) > 1.0? (low_ + up_) / 2.0: est_;
    }
	~MetricEst()
    {}
    double low_, up_, est_, diff_;
};
#endif