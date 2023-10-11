// measure.h
#ifndef MEASURE_H
#define MEASURE_H
class Measure
{
public:
    Measure(): precision_(0.0), recall_(0.0), f1_(0.0), accuracy_(0.0)
    {}
    Measure(const double& precision, const double& recall, const double& f1, const double& NDCG, const double& accuracy): precision_(precision), recall_(recall), f1_(f1), NDCG_(NDCG), accuracy_(accuracy)
    {}
    ~Measure()
    {}
    Measure operator+(const Measure& b)
    {
        Measure measure;
        measure.precision_ = this->precision_ + b.precision_;
        measure.recall_ = this->recall_ + b.recall_;
        measure.f1_ = this->f1_ + b.f1_;
        measure.NDCG_ = this->NDCG_ + b.NDCG_;
        measure.accuracy_ = this->accuracy_ + b.accuracy_;
        return measure;
    }
    Measure & operator+=(const Measure& b)
    {
        this->precision_ += b.precision_;
        this->recall_ += b.recall_;
        this->f1_ += b.f1_;
        this->NDCG_ += b.NDCG_;
        this->accuracy_ += b.accuracy_;
        return *this;
    }
    Measure & operator/=(const int& divisor)
    {
        precision_ = precision_ / divisor;
        recall_ = recall_ / divisor;
        f1_ = f1_ / divisor;
        NDCG_ = NDCG_ / divisor;
        accuracy_ = accuracy_ / divisor;
        return *this;
    }
    double precision_, recall_, f1_, NDCG_, accuracy_;
};
#endif