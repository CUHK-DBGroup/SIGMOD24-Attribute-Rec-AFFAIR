// test.cpp
#include "config.h"
#include "query.h"
#include <chrono>
#include <random>
#include <string>
int main(int argc, char **argv)
{
	Config config;
	for (int i = 0; i < argc; i++) {
        std::string help_str = ""
			"main action [options]\n"
			"\n"
			"options: \n"
			"  --prefix <prefix>\n"
			"  --dataset <dataset>\n"
			"  --support_file <support_file>\n"
			"  --abs_error <error_bound>\n"
			"  --topk <topk>\n"
            "  --initial_size <intial_size>\n"
            "  --query_num <query_num>\n"
            "  --metric <metric_function>\n";
        if (std::string(argv[i]) == "--help") {
            std::cout << help_str << std::endl;
            exit(0);
        }
	}
    config.action_ = argv[1];
    std::cout << "action: " << config.action_ << std::endl;
    for (int i = 0; i < argc; i++) {
        std::string arg = argv[i];
        if (std::string(argv[i]) == "--prefix") {
            config.prefix_ = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--dataset") {
            config.dataset_ = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--predicate1_low") {
            config.predicate1_range_ratio_.first = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--predicate1_up") {
            config.predicate1_range_ratio_.second = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--predicate2_low") {
            config.predicate2_range_ratio_.first = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--predicate2_up") {
            config.predicate2_range_ratio_.second = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--abs_error") {
            config.abs_error_ = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--topk"){
            config.topk_ = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--initial_size"){
            config.initial_size_ = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--batch_size"){
            config.batch_size_ = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--query_num"){
            config.query_num_ = std::stoi(argv[i + 1]);
        }       
        else if (std::string(argv[i]) == "--metric"){
            config.metric_ = argv[i + 1];
            std::cout << config.metric_ << std::endl;
        }   
        else if (arg.substr(0, 2) == "--") {
            std::cerr << "command not recognize " << arg << std::endl;
            exit(1);
        }
    }
    if (config.action_ == "save_serialize_file")
    {
        std::ostringstream oss_data_mat, oss_support_size;
	    oss_data_mat << config.prefix_ << config.dataset_ << ".csv";
        std::vector<std::vector<uint32_t>> data_mat;
        SaveSerializeFile(oss_data_mat.str(), data_mat, config);
    }
    else if (config.action_ == "test")
    {
        printf("test!\n");
        int vector_size = 100000000;
        std::vector<int> input_vector_1(vector_size), input_vector_2(vector_size);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::binomial_distribution<> d1(999, 0.8);
        std::binomial_distribution<> d2(999, 0.75);
        for (int i = 0; i < vector_size; ++i)
        {
            input_vector_1[i] = d1(gen);
            input_vector_2[i] = d2(gen);
        }
        std::vector<std::vector<int>> in_vecs = {input_vector_1, input_vector_2};
        KsTestProcedure(in_vecs, config);
    }
    else if (config.action_ == "topk_query") {
        printf("Topk! Metric--%s\n", config.metric_.c_str());
        std::ostringstream oss_data_mat, oss_support_size;
	    oss_data_mat << config.prefix_ << config.dataset_ << "_predicate_serialize.file";
        oss_support_size << config.prefix_ << config.dataset_ << "_support.csv";
        std::cout << "dataset is " << oss_data_mat.str() << std::endl;
        std::cout << "support is " << oss_support_size.str() << std::endl;
        std::vector<u_int32_t> support_size_vec;
        std::vector<std::vector<uint32_t>> data_mat;
        ReadSupportSize(oss_support_size.str(), support_size_vec);
        LoadSerializeFile(oss_data_mat.str(), data_mat);
        uint32_t num_row = data_mat[0].size();
        printf("supportSizeVec %d dataMat %d %d\n", int(support_size_vec.size()), int(data_mat.size()), num_row);
        config.fail_pr_ = 1.0 / num_row;
        config.predicate1_range_row_ = std::pair<uint32_t, uint32_t>(round(num_row * config.predicate1_range_ratio_.first), round(num_row * config.predicate1_range_ratio_.second));
        config.predicate2_range_row_ = std::pair<uint32_t, uint32_t>(round(num_row * config.predicate2_range_ratio_.first), round(num_row * config.predicate2_range_ratio_.second));
        config.predicate_size_.first = config.predicate1_range_row_.second - config.predicate1_range_row_.first;
        config.predicate_size_.second = config.predicate2_range_row_.second - config.predicate2_range_row_.first;
        printf("predicate1 %f %f %d %d %d\n", config.predicate1_range_ratio_.first, config.predicate1_range_ratio_.second, config.predicate1_range_row_.first, config.predicate1_range_row_.second, config.predicate_size_.first);
        printf("predicate2 %f %f %d %d %d\n", config.predicate2_range_ratio_.first, config.predicate2_range_ratio_.second, config.predicate2_range_row_.first, config.predicate2_range_row_.second, config.predicate_size_.second);
        TopkQuery(data_mat, support_size_vec, config);
    }
    else {
        printf("Command Error!\n");
    }
    return 0;
}