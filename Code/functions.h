// data.h
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <chrono>
#include <ratio>
#include "ks_est.h"
#include "relative_error.h"
#include "absolute_error.h"
#include "pred_pair.h"
#include "config.h"
#include <random>

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

bool cmp_low(const Feature& a, const Feature& b)
{
    return a.low_ > b.low_;
}

bool cmp_ptr_low(const Feature* a, const Feature* b)
{
    return a->low_ > b->low_;
}

bool cmp_up(const Feature& a, const Feature& b)
{
    return a.up_ > b.up_;
}

bool cmp_ptr_up(const Feature* a, const Feature* b)
{
    return a->up_ > b->up_;
}


bool cmp_avg(const Feature& a, const Feature& b)
{
    return a.avg_ > b.avg_;
}

bool cmp_est(const Feature& a, const Feature& b)
{
    return a.est_ > b.est_;
}

bool cmp_alt(const Feature& a, const Feature& b)
{
    return a.alt_ > b.alt_;
}

bool cmp_diff(const Feature& a, const Feature& b)
{
    return a.diff_ < b.diff_;
}

bool cmp_low_up_diff(const std::pair<double, double>& a, const std::pair<double, double>& b)
{
    return (a.second - a.first) < (b.second - b.first);
}

bool cmp_vec(const int a, const int b) 
{
    return a < b;
}

bool cmp_idx_value_pair(const std::pair<int, double>& a, const std::pair<int, double>& b)
{
    return a.first < b.first;
}

void ReadFile(std::string file_name, std::vector<std::vector<uint32_t>>& file_array)
{
    std::vector<std::vector<uint32_t>> temp_file_array;
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    getline(in_file, line_str);
    int count = 1;
    while (getline(in_file, line_str))
    {
        count++;
        std::stringstream ss(line_str);
        std::string str;
        std::vector<uint32_t> line_array;  
        while (getline(ss, str, ','))
            line_array.push_back(stod(str));
        temp_file_array.push_back(line_array);
    }
    in_file.close();
    file_array.resize(temp_file_array[0].size());
    for (int i = 0; i < file_array.size(); ++i)
    {
        file_array[i].resize(temp_file_array.size());
    }
    for (int i = 0; i < file_array.size(); ++i)
    {
        for (int j = 0; j < file_array[i].size(); ++j)
        {
            file_array[i][j] = temp_file_array[j][i];
        } 
    }
}

void ReadSupportSize(std::string file_name, std::vector<uint32_t>& support_size_array)
{
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    while (getline(in_file, line_str))
    {
        support_size_array.push_back(stoi(line_str));
    }
}

void RepermutateMat(uint32_t seed, std::vector<std::vector<uint32_t>>& data_mat)
{
    printf("Random Seed %d\n", seed);
    const uint32_t num_row = data_mat[0].size(); 
    std::default_random_engine e(seed);
    std::vector<uint32_t> random_permutation(num_row);
    for (uint32_t i = 0; i < num_row; i++)
    {
        random_permutation[i] = i;
    }
    std::shuffle(random_permutation.begin(), random_permutation.end(), e);
    std::vector<uint32_t> temp_vector(num_row);
    for (uint32_t i = 0; i < data_mat.size(); i++)
    {
        for (uint32_t j = 0; j < num_row; j++)
        {
            temp_vector[j] = data_mat[i][random_permutation[j]];
        }
        data_mat[i] = temp_vector;
    }
    printf("Permutation ");
    for (uint32_t i = 0; i < 10; i++)
    {
        printf("%d ", random_permutation[i]);
    }
    printf("\n");
}

double AccuracyFilter(const uint32_t& num_pair, std::vector<uint32_t> exact, std::vector<uint32_t> est)
{
    std::vector<uint32_t> diff;
    sort(exact.begin(), exact.end());
    sort(est.begin(), est.end());
    set_symmetric_difference(exact.begin(), exact.end(), est.begin(), est.end(), inserter(diff, diff.end()));
    return 1.0 - double(diff.size()) / double(num_pair);
}

double Precision(std::vector<uint32_t> exact, std::vector<uint32_t> est)
{
    if (est.empty())
    {
        if (exact.empty())
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    std::vector<int> true_pos;
    sort(exact.begin(), exact.end());
    sort(est.begin(), est.end());
    set_intersection(exact.begin(), exact.end(), est.begin(), est.end(), inserter(true_pos, true_pos.end()));
    return double(true_pos.size()) / double(est.size());
}

double Recall(std::vector<uint32_t> exact, std::vector<uint32_t> est)
{
    if (exact.empty())
    {
        return 1.0;   
    }
    std::vector<int> true_pos;
    sort(exact.begin(), exact.end());
    sort(est.begin(), est.end());
    set_intersection(exact.begin(), exact.end(), est.begin(), est.end(), inserter(true_pos, true_pos.end()));
    return double(true_pos.size()) / double(exact.size());
}

double F1Measure(const std::vector<uint32_t>& exact, const std::vector<uint32_t>& est)
{
    double precision = Precision(exact, est), recall = Recall(exact, est);
    if (precision <= 0 && recall <= 0)
    {
        return 0.0;
    }
    else
    {
        return 2.0 * precision * recall / (precision + recall);
    }
}

RelativeError EpsilonDeltaRelativeError(const std::vector<double>& exact_topk_value, const std::vector<double>& approx_topk_value, const double& empirical_threshold)
{
    if (exact_topk_value.size() != approx_topk_value.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_topk_value.size()), int(approx_topk_value.size()));
        return RelativeError(0.0, 0.0);
    }
    double max_relative_error = 0.0, total_relative_error = 0.0;
    for (int i = 0; i < exact_topk_value.size(); i++)
    {
        double relative_error;
        if (exact_topk_value[i] >= empirical_threshold)
        {
            relative_error = abs(exact_topk_value[i] - approx_topk_value[i]) / exact_topk_value[i];
        }
        else
        {
            relative_error = abs(exact_topk_value[i] - approx_topk_value[i]) / empirical_threshold;
        }
        if (relative_error > max_relative_error)
        {
            max_relative_error = relative_error;
        }
        total_relative_error += relative_error;
    }
    return RelativeError(max_relative_error, total_relative_error / exact_topk_value.size());
}

AbsoluteError SelectAbsoluteError(const std::vector<double>& exact_variance, const std::vector<double>& approx_variance)
{
    if (exact_variance.size() != approx_variance.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_variance.size()), int(approx_variance.size()));
        return AbsoluteError(0.0, 0.0);
    }
    double max_absolute_error = 0.0, total_absolute_error = 0.0;
    for (int i = 0; i < exact_variance.size(); i++)
    {
        double absolute_error = abs(exact_variance[i] - approx_variance[i]);
        if (absolute_error > max_absolute_error)
        {
            max_absolute_error = absolute_error;
        }
        total_absolute_error += absolute_error;
    }
    return AbsoluteError(max_absolute_error, total_absolute_error / exact_variance.size());
}

RelativeError TopkRelativeError(const std::vector<double>& exact_topk_value, const std::vector<double>& approx_topk_value)
{
    if (exact_topk_value.size() != approx_topk_value.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_topk_value.size()), int(approx_topk_value.size()));
        return RelativeError(0.0, 0.0);
    }
    double max_relative_error = 0.0, mean_relative_error = 0.0;
    for (int i = 0; i < exact_topk_value.size(); i++)
    {
        if (exact_topk_value[i] > 0.0)
        {
            if (abs(exact_topk_value[i] - approx_topk_value[i]) / exact_topk_value[i] > max_relative_error)
            {
                max_relative_error = abs(exact_topk_value[i] - approx_topk_value[i]) / exact_topk_value[i];
            }
            mean_relative_error += abs(exact_topk_value[i] - approx_topk_value[i]) / exact_topk_value[i];
        }
        else if (exact_topk_value[i] == 0.0)
        {
            if (exact_topk_value[i] - approx_topk_value[i] == 0.0)
            {
                mean_relative_error += 0.0;
            }
            else 
            {
                printf("Try to divide 0!\n");
                printf("i %d value %f\n", i, exact_topk_value[i]);
            }
        }
        else
        {
            printf("Try to divide negative value!\n");
        }
    }
    return RelativeError(max_relative_error, mean_relative_error / exact_topk_value.size());
}

AbsoluteError TopkAbsoluteError(const std::vector<double>& exact_topk_value, const std::vector<double>& approx_topk_value)
{
    if (exact_topk_value.size() != approx_topk_value.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_topk_value.size()), int(approx_topk_value.size()));
        return AbsoluteError(0.0, 0.0);
    }
    double max_absolute_error = 0.0, total_absolute_error = 0.0;
    for (int i = 0; i < exact_topk_value.size(); i++)
    {
        double absolute_error = abs(exact_topk_value[i] - approx_topk_value[i]);
        if (absolute_error > max_absolute_error)
        {
            max_absolute_error = absolute_error;
        }
        total_absolute_error += absolute_error;
    }
    return AbsoluteError(max_absolute_error, total_absolute_error / exact_topk_value.size());
}

RelativeError FilterRelativeError(const std::vector<std::pair<uint32_t, double>>& exact_filter_pair, const std::vector<std::pair<uint32_t, double>>& approx_filter_pair)
{
    if (exact_filter_pair.size() != approx_filter_pair.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_filter_pair.size()), int(approx_filter_pair.size()));
        return RelativeError(0.0, 0.0);
    }
    int cnt = 0;
    double max_relative_error = 0.0, mean_relative_error = 0.0;
    std::vector<std::pair<uint32_t, double>> temp_exact_pair(exact_filter_pair), temp_approx_pair(approx_filter_pair);
    sort(temp_exact_pair.begin(), temp_exact_pair.end(), cmp_idx_value_pair);
    sort(temp_approx_pair.begin(), temp_approx_pair.end(), cmp_idx_value_pair);
    for (int i = 0; i < exact_filter_pair.size(); i++)
    {
        if (temp_exact_pair[i].first == temp_approx_pair[i].first)
        {
            if (abs(temp_exact_pair[i].second - temp_approx_pair[i].second) / temp_exact_pair[i].second > max_relative_error)
            {
                max_relative_error = abs(temp_exact_pair[i].second - temp_approx_pair[i].second) / temp_exact_pair[i].second;
            }
            mean_relative_error += abs(temp_exact_pair[i].second - temp_approx_pair[i].second) / temp_exact_pair[i].second;
            cnt ++;
        }
        else
        {
            continue;
        }
    }
    return RelativeError(max_relative_error, mean_relative_error / cnt);
}

void MetricRank(const std::vector<std::vector<uint32_t>>& data_mat, const std::vector<uint32_t>& support_vec, const Config& config, std::vector<uint32_t>& topk_idx, std::vector<double>& topk_value, double& predicate_time)
{
    uint32_t num_feature = support_vec.size(), num_element = data_mat[0].size(), iteration = 0;
    uint64_t count = 0;

    double pf0 = config.fail_pr_ * 1.0 / (num_feature * (num_element / config.batch_size_ + 1));
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> sample_hist_vec(num_feature);
    std::vector<std::vector<uint32_t>> data_hist(num_feature);

    for (uint32_t i = 0; i < num_feature; ++i)
    {
        data_hist[i].resize(support_vec[i]);
        sample_hist_vec[i].first.resize(support_vec[i], 0);
        sample_hist_vec[i].second.resize(support_vec[i], 0);
    }
    uint32_t num_batch = ceil(1.0 * (num_element - config.initial_size_) / config.batch_size_) + 1;
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> predicate_sample_index_pair_vec(num_batch);
    std::vector<std::pair<uint32_t, uint32_t>> predicate_sample_size_vec(num_batch, std::pair<uint32_t, uint32_t>{0, 0});

    std::vector<Feature> feature_vec, temp(num_feature);
    std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
    predicate_time = 0.0;
    start = std::chrono::high_resolution_clock::now();
    for (uint32_t i = 0; i < config.initial_size_; i++)
    {
        if (data_mat[num_feature][i] >= config.predicate1_range_row_.first && data_mat[num_feature][i] < config.predicate1_range_row_.second)
        {
            predicate_sample_index_pair_vec[0].first.emplace_back(i);
            predicate_sample_size_vec[0].first++;
            // test_count ++;
        } 
        if (data_mat[num_feature][i] >= config.predicate2_range_row_.first && data_mat[num_feature][i] < config.predicate2_range_row_.second)
        {
            predicate_sample_index_pair_vec[0].second.emplace_back(i);
            predicate_sample_size_vec[0].second++;
        } 
    }
    end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    predicate_time += duration.count();
    printf("predicate 1 current size %d cumulative size %d\n", predicate_sample_index_pair_vec[0].first.size(), predicate_sample_size_vec[0].first);
    printf("predicate 2 current size %d cumulative size %d\n", predicate_sample_index_pair_vec[0].second.size(), predicate_sample_size_vec[0].second);
    for (uint32_t i = 0; i < num_feature; i++)
    {
        MetricEst metric("baseline", config, data_mat[i], predicate_sample_index_pair_vec[0], predicate_sample_size_vec[0], support_vec[i], pf0, sample_hist_vec[i], i);
        feature_vec.emplace_back(i, metric.low_, metric.up_, metric.est_, config.initial_size_);
    }

    count += num_feature * config.initial_size_;
    std::set<uint32_t> topk_est, topk_alt;
    do
    {
        std::pair<uint32_t, uint32_t> sample_pair;
        topk_est.clear();
        topk_alt.clear();
        std::set<uint32_t> topk_diff;
        temp.assign(feature_vec.begin(), feature_vec.end());
        std::nth_element(temp.begin(), temp.begin() + config.topk_ - 1, temp.end(), cmp_est);
        sort(temp.begin(), temp.begin() + config.topk_, cmp_est);
        for (int i = 0; i < config.topk_; i++)
        {
            topk_est.insert(temp[i].idx_);
            temp[i].alt_ = temp[i].low_;
        }

        printf("Num of Iteration: %d\n", iteration);
        printf("Est\n");
        for (int i = 0; i < config.topk_ * 2; ++i)
        {
            printf("idx %3d low %2.6f up %1.6f est %1.6f alt %1.6f diff %1.10f sample %d\n", temp[i].idx_, temp[i].low_, temp[i].up_, temp[i].est_, temp[i].alt_, temp[i].up_ - temp[i].low_, temp[i].sample_size_);
        } 

        std::nth_element(temp.begin(), temp.begin() + config.topk_ - 1, temp.end(), cmp_alt);
        sort(temp.begin(), temp.begin() + config.topk_, cmp_alt);

        for (int i = 0; i < config.topk_; i++)
        {
            topk_alt.insert(temp[i].idx_);
        }

        printf("Alt\n");
        for (int i = 0; i < config.topk_ * 2; i++)
        {
            printf("idx %3d low %2.6f up %1.6f est %1.6f alt %1.6f diff %1.10f sample %d\n", temp[i].idx_, temp[i].low_, temp[i].up_, temp[i].est_, temp[i].alt_, temp[i].up_ - temp[i].low_, temp[i].sample_size_);
        } 

        if (topk_est == topk_alt)
        {
            sort(temp.begin(), temp.begin() + config.topk_, cmp_alt);
            topk_idx.clear();
            topk_value.clear();
            printf("Ks Rank has been found.\n");
            for (uint32_t i = 0; i < config.topk_; i++)
            {
                printf("idx %3d low %1.6f up %1.6f est %1.6f diff %1.10f sample %d\n", temp[i].idx_, temp[i].low_, temp[i].up_, temp[i].est_, temp[i].diff_, temp[i].sample_size_);
                topk_idx.push_back(temp[i].idx_);
                topk_value.emplace_back(temp[i].est_);
            }
            printf("PredicateTime %f\n", predicate_time);
            break;
        }
        set_symmetric_difference(topk_est.begin(), topk_est.end(), topk_alt.begin(), topk_alt.end(), inserter(topk_diff, topk_diff.begin()));
        uint32_t idx_max_diff;
        double max_diff = 0;
        for (std::set<uint32_t>::iterator iter = topk_diff.begin(); iter != topk_diff.end(); ++iter)
        {
            if (feature_vec[*iter].diff_ > max_diff)
            {
                max_diff = feature_vec[*iter].up_ - feature_vec[*iter].low_;
                idx_max_diff = *iter;
            }
        }
        if (feature_vec[idx_max_diff].sample_size_ + config.batch_size_ <= num_element)
        {
            iteration ++;
            sample_pair.first = feature_vec[idx_max_diff].sample_size_;
            feature_vec[idx_max_diff].sample_size_ += config.batch_size_;
            sample_pair.second = feature_vec[idx_max_diff].sample_size_;
        }
        else if (max_diff > 0.0)
        {
            iteration ++;
            sample_pair.first = feature_vec[idx_max_diff].sample_size_;
            feature_vec[idx_max_diff].sample_size_ = num_element;
            sample_pair.second = feature_vec[idx_max_diff].sample_size_;
        }
        else
        {
            std::cout << "Elements have been used up!" << std::endl;
            std::nth_element(temp.begin(), temp.begin() + config.topk_ - 1, temp.end(), cmp_est);
            sort(temp.begin(), temp.begin() + config.topk_, cmp_est);
            topk_idx.clear();
            topk_value.clear();
            printf("Ks Rank has been found.\n");
            for (uint32_t i = 0; i < config.topk_; ++i)
            {
                printf("idx %3d low %1.6f up %1.6f est %1.6f alt %1.6f diff %1.10f sample %d\n", temp[i].idx_, temp[i].low_, temp[i].up_, temp[i].est_, temp[i].alt_, temp[i].up_ - temp[i].low_, temp[i].sample_size_);
                topk_idx.push_back(temp[i].idx_);
                topk_value.push_back(temp[i].est_);
            }
            break;
        }
        // count += sample_size - prev_sample_size;
        uint32_t current_batch = ceil(1.0 * (sample_pair.second - config.initial_size_) / config.batch_size_);
        if (predicate_sample_index_pair_vec[current_batch].first.empty())
        {
            // important
            start = std::chrono::high_resolution_clock::now();
            predicate_sample_size_vec[current_batch] = predicate_sample_size_vec[current_batch - 1];
            for (uint32_t i = sample_pair.first; i < sample_pair.second; i++)
            {
                if (data_mat[num_feature][i] >= config.predicate1_range_row_.first && data_mat[num_feature][i] < config.predicate1_range_row_.second)
                {
                    predicate_sample_index_pair_vec[current_batch].first.emplace_back(i);
                    predicate_sample_size_vec[current_batch].first++;
                } 
                if (data_mat[num_feature][i] >= config.predicate2_range_row_.first && data_mat[num_feature][i] < config.predicate2_range_row_.second)
                {
                    predicate_sample_index_pair_vec[current_batch].second.emplace_back(i);
                    predicate_sample_size_vec[current_batch].second++;
                } 
            }
            end = std::chrono::high_resolution_clock::now();
	        duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            predicate_time += duration.count();
            // std::sort(predicate_sample_index_pair_vec[current_batch].first.begin(), predicate_sample_index_pair_vec[current_batch].first.end());
            // std::sort(predicate_sample_index_pair_vec[current_batch].second.begin(), predicate_sample_index_pair_vec[current_batch].second.end());
        }
        printf("sample size %d current batch %d\n", sample_pair.second, current_batch);
        printf("predicate 1 current size %d cumulative size %d\n", predicate_sample_index_pair_vec[current_batch].first.size(), predicate_sample_size_vec[current_batch].first);
        printf("predicate 2 current size %d cumulative size %d\n", predicate_sample_index_pair_vec[current_batch].second.size(), predicate_sample_size_vec[current_batch].second);

        MetricEst metric("baseline", config, data_mat[idx_max_diff], predicate_sample_index_pair_vec[current_batch], predicate_sample_size_vec[current_batch], support_vec[idx_max_diff], pf0, sample_hist_vec[idx_max_diff], idx_max_diff);
        feature_vec[idx_max_diff].UpdateWithSampleSize(metric.low_, metric.up_, metric.diff_, metric.est_, sample_pair.second);
    } while (topk_est != topk_alt);
}

double ExactCd(const std::vector<uint32_t>& in_vec1, const std::vector<uint32_t>& in_vec2, const int support_size) 
{
    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> in_hist_pair(std::vector<uint32_t>(support_size, 0), std::vector<uint32_t>(support_size, 0));
    
    auto update_hist = [](const std::vector<uint32_t>& in_vec, std::vector<uint32_t>& in_hist) -> void 
    {
        for (uint32_t i = 0; i < in_vec.size(); i++) {
            in_hist[in_vec[i]] ++;
        }
    };
    update_hist(in_vec1, in_hist_pair.first);
    update_hist(in_vec2, in_hist_pair.second);
    double max_diff_prob = 0.0;
    for (uint32_t i = 0; i < support_size; ++i) 
    {
        double curr_diff_prob = abs(1.0 * in_hist_pair.first[i] / in_vec1.size() - 1.0 * in_hist_pair.second[i] / in_vec2.size());
        if (curr_diff_prob > max_diff_prob) {
            max_diff_prob = curr_diff_prob;
        }
    }
    return max_diff_prob;
}

double ExactKs(const std::vector<uint32_t>& in_vec1, const std::vector<uint32_t>& in_vec2, const int support_size) 
{
    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> in_hist_pair(std::vector<uint32_t>(support_size, 0), std::vector<uint32_t>(support_size, 0));
    auto update_hist = [](const std::vector<uint32_t>& in_vec, std::vector<uint32_t>& in_hist) -> void 
    {
        for (uint32_t i = 0; i < in_vec.size(); i++) 
        {
            in_hist[in_vec[i]] ++;
        }
    };
    update_hist(in_vec1, in_hist_pair.first);
    update_hist(in_vec2, in_hist_pair.second);
    double max_diff_prob = 0.0;
    std::pair<uint32_t, uint32_t> cum_hist_pair(0, 0);
    for (uint32_t i = 0; i < support_size; ++i) 
    {
        cum_hist_pair.first += in_hist_pair.first[i];
        cum_hist_pair.second += in_hist_pair.second[i];
        double curr_diff_prob = abs(1.0 * cum_hist_pair.first / in_vec1.size() - 1.0 * cum_hist_pair.second / in_vec2.size());
        if (curr_diff_prob > max_diff_prob) 
        {
            max_diff_prob = curr_diff_prob;
        }
    }
    return max_diff_prob;
}

double ExactOneNorm(const std::vector<std::vector<int>>& in_vecs, const int support_size) {
    std::vector<std::vector<int>> in_hists(2, std::vector<int>(support_size, 0));
    for (int i = 0; i < in_hists.size(); i++) 
    {
        for (int j = 0; j < in_vecs[i].size(); j++)
        {
            in_hists[i][in_vecs[i][j]] ++;
        }
    }
    double max_diff_prob = 0.0;
    for (int i = 0; i < support_size; ++i) 
    {
        double curr_diff_prob = abs(1.0 * in_hists[0][i] / in_vecs[0].size() - 1.0 * in_hists[1][i] / in_vecs[1].size());
        max_diff_prob += curr_diff_prob;
    }
    return max_diff_prob;
}

void ExpectedBoundIntegral(const std::vector<std::vector<int>>& in_vecs, std::vector<std::vector<int>>& sample_hists, const int sample_size, const int support_size, std::vector<double>& expected_diff_pair)
{
    const double init_fair_prob = 0.5;
    for (int i = 0; i < 2; i++)
    {
        double expected_diff = 0.0, single_fail_prob = init_fair_prob;
        while (single_fail_prob >= 1.0 / in_vecs[i].size())
        {
            double max_diff = 1.0;
            double a = log(2.0 / single_fail_prob);
            for (int j = 0; j < sample_hists[i].size(); j++)
            {
                double est_x = 1.0 * sample_hists[i][j] / sample_size, 
                lower_bound = std::max(0.0, (est_x + 2.0 * a / 3.0 / sample_size - sqrt(4.0 * a * a / 3.0 / sample_size / sample_size * (1.0 / 3.0 + est_x) + 2.0 * a / sample_size * est_x * (1.0 - est_x))) / (2.0 * a / sample_size + 1.0)),
                upper_bound = std::min(1.0, est_x + 1.0 * a / sample_size + sqrt(2.0 * a / sample_size * (est_x + a / 2.0 / sample_size))),
                curr_diff = upper_bound - lower_bound;
                if (curr_diff < max_diff)
                {
                    max_diff = curr_diff;
                }
            }
            expected_diff += max_diff * single_fail_prob;
            single_fail_prob /= 2;
        }
        expected_diff += 1.0 * single_fail_prob;
        expected_diff_pair[i] = expected_diff; 
    }
}

void ApproxKsTest(const std::vector<std::vector<int>>& in_vecs, const int support_size, const double fail_prob) 
{
    int max_sample_size = 100000000;
    std::pair<int, int> sample_pair(0, 1024);
    std::vector<std::vector<int>> sample_hists(2, std::vector<int>(support_size, 0));
    std::vector<int> in_vec_size = {int(in_vecs[0].size()), int(in_vecs[1].size())};
    while (sample_pair.second <= max_sample_size)
    {
        double rel_diff_integral = 0.0, rel_diff_baseline = 0.0;
        for (int i = 0; i < 2; i++)
        {
            for (int j = sample_pair.first; j < sample_pair.second; j++)
            {
                sample_hists[i][in_vecs[i][j]] ++;
            }  
        }
        double est_ks_test = 0.0;
        for (int i = 0; i < support_size; i++)
        {
            double curr_diff_prob = abs(1.0 * sample_hists[0][i] / sample_pair.second - 1.0 * sample_hists[1][i] / sample_pair.second);
            // abs_diff_vec[i] = curr_diff_prob;
            if (est_ks_test < curr_diff_prob)
            {
                est_ks_test = curr_diff_prob;
            }         
        }
        printf("current sample_size is %d est_ks_test is %f\n", sample_pair.second, est_ks_test);
        std::vector<double> expected_bound(2, 0.0);
        ExpectedBoundIntegral(in_vecs, sample_hists, sample_pair.second, support_size, expected_bound);
        for (int i = 0; i < 2; i++) 
        {
            rel_diff_integral += expected_bound[i] + sqrt(1.0 * (in_vec_size[i] - sample_pair.second) / (1.0 * sample_pair.second * in_vec_size[i])) * sqrt(2.0 * log(2.0 / fail_prob));
            rel_diff_baseline += 1.0 / 2.0 * sqrt((support_size - 1.0) * ((in_vec_size[i] - sample_pair.second) / (1.0 * sample_pair.second * (in_vec_size[i] - 1.0)))) + sqrt(1.0 * (in_vec_size[i] - sample_pair.second) / (1.0 * sample_pair.second * in_vec_size[i])) * sqrt(2.0 * log(2.0 / fail_prob));
            printf("test result is %f %f %f\n", 1.0 / 2.0 * sqrt((support_size - 1.0) * ((in_vec_size[i] - sample_pair.second) / (1.0 * sample_pair.second * (in_vec_size[i] - 1.0)))), sqrt(1.0 * (in_vec_size[i] - sample_pair.second) / (1.0 * sample_pair.second * in_vec_size[i])) * sqrt(2.0 * log(2.0 / fail_prob)), 1.0 / 2.0 * sqrt((support_size - 1.0) * ((in_vec_size[i] - sample_pair.second) / (1.0 * sample_pair.second * (in_vec_size[i] - 1.0)))) + sqrt(1.0 * (in_vec_size[i] - sample_pair.second) / (1.0 * sample_pair.second * in_vec_size[i])) * sqrt(2.0 * log(2.0 / fail_prob)));
        }
        printf("integral bound is %f\n", rel_diff_integral);
        printf("baseline bound is %f\n", rel_diff_baseline);
        printf("integral sample_size %d low %f up %f\n", sample_pair.second, std::max(0.0, est_ks_test - rel_diff_integral), std::min(1.0, est_ks_test + rel_diff_integral));    
        printf("baseline sample_size %d low %f up %f\n", sample_pair.second, std::max(0.0, est_ks_test - rel_diff_baseline), std::min(1.0, est_ks_test + rel_diff_baseline));  
        sample_pair.first = sample_pair.second;
        if (sample_pair.second > max_sample_size / 2 && sample_pair.second < max_sample_size) {
            sample_pair.second = max_sample_size;
        }
        else {
            sample_pair.second *= 2;
        }
    }
}

double ExactMetric(const Config& config, const std::vector<uint32_t>& in_vec, const std::pair<std::vector<uint32_t>, std::vector<uint32_t>>& predicate_sample_index_pair, const int support_size)
{
    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> in_hist_pair(std::vector<uint32_t>(support_size, 0), std::vector<uint32_t>(support_size, 0));
    auto update_hist = [](const std::vector<uint32_t>& in_vec, std::vector<uint32_t>& in_hist, const std::vector<uint32_t>& predicate_sample_index) -> void
    {
        for (uint32_t i = 0; i < predicate_sample_index.size(); i++) 
        {
            in_hist[in_vec[predicate_sample_index[i]]] ++;
        }  
    };
    update_hist(in_vec, in_hist_pair.first, predicate_sample_index_pair.first);
    update_hist(in_vec, in_hist_pair.second, predicate_sample_index_pair.second);
    double returned_metric_function= 0.0;

    if (config.metric_ == "ks") 
    {
        std::pair<uint32_t, uint32_t> cum_hist_pair(0, 0);
        for (uint32_t i = 0; i < support_size; ++i) 
        {
            cum_hist_pair.first += in_hist_pair.first[i];
            cum_hist_pair.second += in_hist_pair.second[i];
            double curr_diff_prob = abs(1.0 * cum_hist_pair.first / config.predicate_size_.first - 1.0 * cum_hist_pair.second / config.predicate_size_.second);
            if (curr_diff_prob > returned_metric_function) 
            {
                returned_metric_function = curr_diff_prob;
            }
        }
    }
    else if (config.metric_ == "cd") 
    {
        for (uint32_t i = 0; i < support_size; ++i) 
        {
            double curr_diff_prob = abs(1.0 * in_hist_pair.first[i] / config.predicate_size_.first - 1.0 * in_hist_pair.second[i] / config.predicate_size_.second);
            if (curr_diff_prob > returned_metric_function) {
                returned_metric_function = curr_diff_prob;
            }
        }
    }
    else if (config.metric_ == "emd") 
    {
        for (uint32_t i = 0; i < support_size; ++i) 
        {
            double curr_diff_prob = abs(1.0 * in_hist_pair.first[i] / config.predicate_size_.first - 1.0 * in_hist_pair.second[i] / config.predicate_size_.second);
            returned_metric_function += curr_diff_prob;
        }
        returned_metric_function /= 2.0;
    }
    else if (config.metric_ == "euclidean") 
    {
        for (uint32_t i = 0; i < support_size; ++i) 
        {
            double curr_diff_prob = abs(1.0 * in_hist_pair.first[i] / config.predicate_size_.first - 1.0 * in_hist_pair.second[i] / config.predicate_size_.second);
            returned_metric_function += curr_diff_prob * curr_diff_prob;
        }
        returned_metric_function = sqrt(returned_metric_function / 2.0);
    }
    else
    {
        printf("unrecognize metric function\n");
    }
    return returned_metric_function;
}

void ExactMetricTopk(const std::vector<std::vector<uint32_t>>& data_mat, const std::vector<uint32_t>& support_vec, std::vector<uint32_t>& topk_idx, std::vector<uint32_t>& topk_idx_plus, std::vector<double>& all_value, std::vector<double>& topk_value, const Config& config, double& predicate_time)
{
	std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();

    uint32_t num_feature = support_vec.size(), num_element = data_mat[0].size();
    if (num_feature != data_mat.size() - 1)
    {
        printf("num_feature has error!\n");
    }
    std::vector<Feature> feature_vec;

    std::chrono::high_resolution_clock::time_point cal_start, cal_end;
    std::chrono::duration<double> cal_duration;
    cal_start = std::chrono::high_resolution_clock::now();


    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> predicate_sample_index_pair;
    predicate_time = 0.0;
    start = std::chrono::high_resolution_clock::now();
    for (uint32_t i = 0; i < num_element; i++)
    {
        if (data_mat[num_feature][i] >= config.predicate1_range_row_.first && data_mat[num_feature][i] < config.predicate1_range_row_.second)
        {
            predicate_sample_index_pair.first.emplace_back(i);
        } 
        if (data_mat[num_feature][i] >= config.predicate2_range_row_.first && data_mat[num_feature][i] < config.predicate2_range_row_.second)
        {
            predicate_sample_index_pair.second.emplace_back(i);
        } 
    }
    end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    predicate_time += duration.count();
    for (uint32_t i = 0; i < num_feature; i++)
    {
        double exact_variance_score = ExactMetric(config, data_mat[i], predicate_sample_index_pair, support_vec[i]);
        feature_vec.emplace_back(i, exact_variance_score);
    }
    cal_end = std::chrono::high_resolution_clock::now();
	cal_duration = std::chrono::duration_cast<std::chrono::duration<double>>(cal_end - cal_start);
    printf("CalDuration %f\n", cal_duration.count());

    topk_idx.clear();
    topk_idx_plus.clear();
    topk_value.clear();
    std::vector<Feature> temp(feature_vec);
    sort(temp.begin(), temp.end(), cmp_est);
    for (int i = 0; i < temp.size(); i++)
    {
        if (i < config.topk_) {
            topk_idx.push_back(temp[i].idx_);
            topk_idx_plus.push_back(temp[i].idx_);  
            topk_value.push_back(temp[i].est_);
        }
        else if (abs(temp[i].est_ - temp[config.topk_ - 1].est_) <= 1e-12) 
        {// With difference smaller than 1e-12, we regard them as equal values.  
            topk_idx_plus.push_back(temp[i].idx_);
        }
        if (i < 20) 
        {
            printf("idx %3d avg %1.6f\n", temp[i].idx_, temp[i].est_);
        }
    }
    all_value.resize(num_feature);
    for (int i = 0; i < feature_vec.size(); i++) 
    {
        all_value[i] = feature_vec[i].est_;
    }
    end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    printf("Feature %d Count %ld ExactDuration %f PredicateTime %f\n", num_feature, long(num_feature) * long(data_mat[0].size()), duration.count(), predicate_time);
}

void ApproxMetricTopk(const std::vector<std::vector<uint32_t>>& data_mat, const std::vector<uint32_t>& support_vec, const Config& config, std::vector<uint32_t>& topk_idx, std::vector<double>& topk_value, double& predicate_time)
{
    const uint32_t num_feature = support_vec.size(), num_element = data_mat[0].size();
    std::pair<uint32_t, uint32_t> sample_pair(0, config.initial_size_);

    const double pf0 = config.fail_pr_ * 1.0 / ceil(log2(num_element * 1.0 / config.initial_size_)); //require_diff = config.abs_error_ * 2
    double access_data_time = 0.0, sample_index_time = 0.0;
    printf("Element %d Feature %d pf0 %1.15f\n", num_element, num_feature, pf0);
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> sample_hist_vec(num_feature);
    std::vector<std::pair<uint32_t, uint32_t>> predicate_sample_size_vec;
    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> predicate_sample_index_pair;

    std::vector<std::pair<Feature*, Feature*>> low_up_vec(config.topk_, std::pair<Feature*, Feature*>(NULL, NULL));
    for (uint32_t i = 0; i < num_feature; i++)
    {
        sample_hist_vec[i].first.resize(support_vec[i], 0);
        sample_hist_vec[i].second.resize(support_vec[i], 0);       
    }
    
    std::vector<Feature*> feature_vec, result_vec;
    
    for (uint32_t i = 0; i < num_feature; i++) 
    {
        feature_vec.emplace_back(new Feature(i, 0.0));
    }
    topk_idx.clear();
    topk_value.clear();
    uint32_t num_iter = 0;

    std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
    predicate_time = 0.0;
    while (!feature_vec.empty() && sample_pair.second <= num_element)
    {
        uint32_t exited_feature = 0; //num_scan = 0, 
        start = std::chrono::high_resolution_clock::now();
        for (uint32_t i = sample_pair.first; i < sample_pair.second; i++)
        {
            if (data_mat[num_feature][i] >= config.predicate1_range_row_.first && data_mat[num_feature][i] < config.predicate1_range_row_.second)
            {
                predicate_sample_index_pair.first.emplace_back(i);
            } 
            if (data_mat[num_feature][i] >= config.predicate2_range_row_.first && data_mat[num_feature][i] < config.predicate2_range_row_.second)
            {
                predicate_sample_index_pair.second.emplace_back(i);
            } 
        }
        
        predicate_sample_size_vec.emplace_back(predicate_sample_index_pair.first.size(), predicate_sample_index_pair.second.size());
        printf("num_iter %d\n", num_iter);
        printf("predicate 1 current size %d cumulative size %d\n", predicate_sample_index_pair.first.size(), predicate_sample_size_vec[num_iter].first);
        printf("predicate 2 current size %d cumulative size %d\n", predicate_sample_index_pair.second.size(), predicate_sample_size_vec[num_iter].second);
        
        end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        predicate_time += duration.count();
    
        for (uint32_t i = 0; i < feature_vec.size(); i++)
        {
            if (feature_vec[i]->diff_ > config.abs_error_ || feature_vec[i]->to_update_ == true)
            {
                exited_feature ++;
                uint32_t feature_idx = feature_vec[i]->idx_;
                uint32_t prev_iter = ceil(log2(feature_vec[i]->sample_size_ / config.initial_size_)), temp_iter = ceil(log2(sample_pair.second / config.initial_size_));
                std::pair<std::pair<uint32_t, uint32_t>, std::pair<uint32_t, uint32_t>> 
                 single_sample_range_pair(std::pair<uint32_t, uint32_t>(predicate_sample_size_vec[prev_iter].first, predicate_sample_size_vec[num_iter].first), 
                 std::pair<uint32_t, uint32_t>(predicate_sample_size_vec[prev_iter].second, predicate_sample_size_vec[num_iter].second));
                MetricEst metric("approx", config, data_mat[feature_idx], predicate_sample_index_pair, single_sample_range_pair, support_vec[feature_idx], pf0, sample_hist_vec[feature_idx], feature_idx);
                feature_vec[i]->UpdateWithSampleSize(metric.low_, metric.up_, metric.diff_, metric.est_, sample_pair.second);
            }
        }


        printf("SampleSize %d ExistedFeature %d\n", sample_pair.second, exited_feature);
        std::nth_element(feature_vec.begin(), feature_vec.begin() + config.topk_ - 1, feature_vec.end(), [](auto a, auto b){return a->low_ > b->low_;});
        double kth_low = feature_vec[config.topk_ - 1]->low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + config.topk_ - 1, feature_vec.end(), [](auto a, auto b){return a->est_ > b->est_;});
        sort(feature_vec.begin(), feature_vec.begin() + config.topk_, [](auto a, auto b){return a->est_ > b->est_;});
        double max_diff = (*std::max_element(feature_vec.begin(), feature_vec.begin() + config.topk_, [](Feature* a, Feature* b){return a->diff_ < b->diff_;}))->diff_;
        for (int i = 0; i < config.topk_; i++)
        {
            printf("idx %3d low %1.6f up %1.6f est %1.6f avg %1.6f diff %1.10f supportSize %d b_low %1.6f b_up %1.6f n_low %1.6f n_up %1.6f\n",
             feature_vec[i]->idx_, feature_vec[i]->low_, feature_vec[i]->up_, feature_vec[i]->est_, feature_vec[i]->avg_, feature_vec[i]->diff_, support_vec[feature_vec[i]->idx_], feature_vec[i]->base_low_, feature_vec[i]->base_up_, feature_vec[i]->new_low_, feature_vec[i]->new_up_);
        } 
        printf("max_diff %f\n", max_diff);
        if (max_diff <= config.abs_error_)
        {
            if (sample_pair.second == num_element)
            {
                printf("ApproximateVarianceTopk has been found when data are used up.\n");
                for (uint32_t i = 0; i < config.topk_; i++)
                {
                    topk_idx.push_back(feature_vec[i]->idx_);
                    topk_value.push_back(feature_vec[i]->est_);
                    printf("idx %3d low %2.6f up %1.6f est %1.6f avg %1.6f diff %1.10f\n", feature_vec[i]->idx_, feature_vec[i]->low_, feature_vec[i]->up_, feature_vec[i]->est_, feature_vec[i]->avg_, feature_vec[i]->diff_);
                }
                printf("SampleSize %d SampleDataTime %f AccessDataTime %f PredicateTime %f\n", sample_pair.second, sample_index_time, access_data_time, predicate_time);
                break;
            }
            else 
            {
                result_vec.assign(feature_vec.begin(), feature_vec.begin() + config.topk_);
                std::nth_element(feature_vec.begin(), feature_vec.begin() + config.topk_ - 1, feature_vec.end(), cmp_ptr_low);
                sort(feature_vec.begin(), feature_vec.begin() + config.topk_, cmp_ptr_low);
                printf("Lower Bounds\n");
                for (uint32_t i = 0; i < config.topk_; i++)
                {
                    low_up_vec[i].first = feature_vec[i];
                    printf("idx = %d bound = %f\n", feature_vec[i]->idx_, feature_vec[i]->low_);
                }

                std::nth_element(feature_vec.begin(), feature_vec.begin() + config.topk_ - 1, feature_vec.end(), cmp_ptr_up);
                sort(feature_vec.begin(), feature_vec.begin() + config.topk_, cmp_ptr_up);
                printf("Upper Bounds\n");
                for (uint32_t i = 0; i < config.topk_; i++)
                {
                    low_up_vec[i].second = feature_vec[i];
                    printf("idx = %d bound = %f\n", feature_vec[i]->idx_, feature_vec[i]->up_);
                }

                bool low_up_flag = true;
                printf("Test Process\n");
                for (uint32_t i = 0; i < config.topk_; i++)
                {
                    printf("i = %d est-low %f up-est %f est %f low %f up %f\n", i, result_vec[i]->est_ - low_up_vec[i].first->low_, low_up_vec[i].second->up_ - result_vec[i]->est_, result_vec[i]->est_, low_up_vec[i].first->low_, low_up_vec[i].second->up_);
                    if (std::max(result_vec[i]->est_ - low_up_vec[i].first->low_, low_up_vec[i].second->up_ - result_vec[i]->est_) > config.abs_error_)
                    {
                        low_up_flag = false;
                        result_vec[i]->to_update_ = true;
                        low_up_vec[i].first->to_update_ = true;
                        low_up_vec[i].second->to_update_ = true;
                        printf("set false i = %d\n", i);
                    }
                }
                if (low_up_flag == true)
                {
                    printf("ApproximateVarianceTopk has been found.\n");
                    for (uint32_t i = 0; i < config.topk_; i++)
                    {
                        topk_idx.push_back(result_vec[i]->idx_);
                        topk_value.push_back(result_vec[i]->est_);
                        printf("idx %3d low %2.6f up %1.6f est %1.6f avg %1.6f diff %1.10f\n", result_vec[i]->idx_, result_vec[i]->low_, result_vec[i]->up_, result_vec[i]->est_, result_vec[i]->avg_, result_vec[i]->diff_);
                    }
                    printf("SampleSize %d SampleDataTime %f AccessDataTime %f PredicateTime %f\n", sample_pair.second, sample_index_time, access_data_time, predicate_time);
                    break;
                }
                else if (sample_pair.second > num_element / 2)
                {
                    sample_pair.first = sample_pair.second;
                    sample_pair.second = num_element; 
                }
                else 
                {
                    sample_pair.first = sample_pair.second;
                    sample_pair.second *= 2; 
                }
            }
        }
        else if (sample_pair.second >= num_element)
        {
            printf("ApproximateVarianceTopk has been found when data are used up.\n");
            for (uint32_t i = 0; i < config.topk_; i++)
            {
                topk_idx.push_back(feature_vec[i]->idx_);
                topk_value.push_back(feature_vec[i]->est_);
                printf("idx %3d low %2.6f up %1.6f est %1.6f avg %1.6f diff %1.10f\n", feature_vec[i]->idx_, feature_vec[i]->low_, feature_vec[i]->up_, feature_vec[i]->est_, feature_vec[i]->avg_, feature_vec[i]->diff_);
            }
            printf("SampleSize %d SampleDataTime %f AccessDataTime %f PredicateTime %f\n", sample_pair.second, sample_index_time, access_data_time, predicate_time);
            break;
        }
        else if (sample_pair.second > num_element / 2)
        {
            sample_pair.first = sample_pair.second;
            sample_pair.second = num_element; 
        }
        else
        {
            sample_pair.first = sample_pair.second;
            sample_pair.second *= 2; 
        }
        std::vector<Feature*>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if ((*iter)->up_ < kth_low) 
            {
                delete *iter;
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
        num_iter++;
    } 
}