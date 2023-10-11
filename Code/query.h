#ifndef QUERY_H
#define QUERY_H
#include "functions.h"
#include "config.h"
#include "measure.h"
#include "file_ctrl.h"
#include <random>
void SaveSerializeFile(const std::string in_file_name, std::vector<std::vector<uint32_t>>& data_mat, const Config& config)
{
    ReadFile(in_file_name, data_mat);
	std::string file_int = config.prefix_ + config.dataset_ + "_serialize.file";
    save_serialized_graph(file_int, data_mat);
}

void LoadSerializeFile(const std::string in_file_name, std::vector<std::vector<uint32_t>>& data_mat)
{
    // std::string file_int = config.prefix_ + config.dataset_ + "_serialize.file";
	load_serialized_graph(in_file_name, data_mat);
}

void KsTestProcedure(const std::vector<std::vector<int>>& in_vecs, const Config& conf) 
{
    ApproxKsTest(in_vecs, 1000, 1.0 / 100000000);
}

void TopkQuery(std::vector<std::vector<uint32_t>>& data_mat, const std::vector<uint32_t>& support_vec, const Config& config)
{
    printf("num_feature_plus_index %d num_record %d\n", uint32_t(support_vec.size() + 1), uint32_t(data_mat[0].size()));
	double exact_time = 0.0, baseline_time = 0.0, approx_time = 0.0, predicate_time = 0.0;
	RelativeError approx_rel_err(0.0, 0.0), baseline_rel_err(0.0, 0.0), tmp_rel_err;  
	Measure approx_measure(0.0, 0.0, 0.0, 0.0, 0.0), baseline_measure(0.0, 0.0, 0.0, 0.0, 0.0), temp_measure;

	printf("Initial Sample Size %d\nTop_%d\n", config.initial_size_, config.topk_);
	std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
    for (uint32_t t = 0; t < config.query_num_; t++)
    {
		printf("\nIteration %d\n", t);
		RepermutateMat(t, data_mat);
		
		std::vector<double> exact_topk_value, exact_all_value, approx_topk_value, baseline_topk_value;
		std::vector<uint32_t> exact_topk_idx, exact_topk_idx_plus, approx_topk_idx, baseline_topk_idx;
		start = std::chrono::high_resolution_clock::now();
		ExactMetricTopk(data_mat, support_vec, exact_topk_idx, exact_topk_idx_plus, exact_all_value, exact_topk_value, config, predicate_time);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		printf("ExactTopkTime: %f\n", duration.count());
		exact_time += (duration.count() - predicate_time);

		start = std::chrono::high_resolution_clock::now();
		ApproxMetricTopk(data_mat, support_vec, config, approx_topk_idx, approx_topk_value, predicate_time);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		temp_measure = Measure(Precision(exact_topk_idx_plus, approx_topk_idx), Recall(exact_topk_idx_plus, approx_topk_idx), F1Measure(exact_topk_idx_plus, approx_topk_idx), 0.0, 0.0);
		printf("Approx %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		printf("ApproxTopkTime: %f\n", duration.count());
		approx_time += (duration.count() - predicate_time);
		approx_measure += temp_measure; 

		if (config.metric_ != "ks" && config.metric_ != "cd") 
		{
			start = std::chrono::high_resolution_clock::now();
			MetricRank(data_mat, support_vec, config, baseline_topk_idx, baseline_topk_value, predicate_time);
			end = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
			temp_measure = Measure(Precision(exact_topk_idx_plus, baseline_topk_idx), Recall(exact_topk_idx_plus, baseline_topk_idx), F1Measure(exact_topk_idx_plus, baseline_topk_idx), 0.0, 0.0);
			printf("Baseline %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
			printf("RankTopkTime: %f\n", duration.count());
			baseline_time += (duration.count() - predicate_time);
			tmp_rel_err = TopkRelativeError(exact_topk_value, baseline_topk_value);
			baseline_rel_err += tmp_rel_err;
			baseline_measure += temp_measure; 
		}
	}
	exact_time /= config.query_num_;
	printf("Exact %f\n", exact_time); 

	approx_time /= config.query_num_;
	approx_rel_err /= config.query_num_;
	approx_measure /= config.query_num_;

	baseline_time /= config.query_num_;
	baseline_rel_err /= config.query_num_;
	baseline_measure /= config.query_num_;

	printf("Approx %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", approx_time, approx_rel_err.max_, approx_rel_err.avg_, approx_measure.precision_, approx_measure.recall_, approx_measure.f1_, approx_measure.NDCG_);
	printf("Baseline %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", baseline_time, baseline_rel_err.max_, baseline_rel_err.avg_, baseline_measure.precision_, baseline_measure.recall_, baseline_measure.f1_, baseline_measure.NDCG_);
}

#endif