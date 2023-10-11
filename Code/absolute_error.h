// absolute_error.h
#ifndef ABSOLUTE_ERROR_H
#define ABSOLUTE_ERROR_H
class AbsoluteError
{
public:
    AbsoluteError(): max_(0.0), avg_(0.0)
    {}
    AbsoluteError(const double& MAE, const double& AAE): max_(MAE), avg_(AAE)
    {}
    ~AbsoluteError()
    {}
    AbsoluteError operator+(const AbsoluteError& b)
    {
        AbsoluteError absolute_error;
        absolute_error.max_ = this->max_ + b.max_;
        absolute_error.avg_ = this->avg_ + b.avg_;
        return absolute_error;
    }
    AbsoluteError & operator+=(const AbsoluteError& b)
    {
        this->max_ += b.max_;
        this->avg_ += b.avg_;
        return *this;
    }
    AbsoluteError & operator/=(const int& divisor)
    {
        this->max_ /= divisor;
        this->avg_ /= divisor;
        return *this;
    }
    double max_, avg_;
};
#endif