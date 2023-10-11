// relative_error.h
#ifndef RELATIVE_ERROR_H
#define RELATIVE_ERROR_H
class RelativeError
{
public:
    RelativeError(): max_(0.0), avg_(0.0)
    {}
    RelativeError(const double& MRE, const double& ARE): max_(MRE), avg_(ARE)
    {}
    ~RelativeError()
    {}
    RelativeError operator+(const RelativeError& b)
    {
        RelativeError relative_error;
        relative_error.max_ = this->max_ + b.max_;
        relative_error.avg_ = this->avg_ + b.avg_;
        return relative_error;
    }
    RelativeError & operator+=(const RelativeError& b)
    {
        this->max_ += b.max_;
        this->avg_ += b.avg_;
        return *this;
    }
    RelativeError & operator/=(const int& divisor)
    {
        this->max_ /= divisor;
        this->avg_ /= divisor;
        return *this;
    }
    double max_, avg_;
};
#endif