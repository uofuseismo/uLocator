#include "uLocator/topography/constant.hpp"

using namespace ULocator::Topography;

class Constant::ConstantImpl
{
public:
    double mElevation{0};
    double mMinimumElevation{0};
    double mMaximumElevation{0};
    bool mHaveTopography{false};
};

/// Constructor
Constant::Constant() :
    pImpl(std::make_unique<ConstantImpl> ())
{
}

/// Copy constructor
Constant::Constant(const Constant &topography)
{
    *this = topography;
}

/// Move constructor
Constant::Constant(Constant &&topography) noexcept
{
    *this = std::move(topography);
}

/// Destructor
Constant::~Constant() = default;

/// Copy assignment
Constant& Constant::operator=(const Constant &topography)
{
    if (&topography == this){return *this;}
    pImpl = std::make_unique<ConstantImpl> (*topography.pImpl);
    return *this;
}

/// Move assignment
Constant& Constant::operator=(Constant &&topography) noexcept
{
    if (&topography == this){return *this;}
    pImpl = std::move(topography.pImpl);
    return *this;
}

/// Reset class and release memory
void Constant::clear() noexcept
{
    pImpl = std::make_unique<ConstantImpl> (); 
}

/// Sets the topography
void Constant::set(const double elevation) noexcept
{
    pImpl->mElevation = elevation;
    pImpl->mMinimumElevation = elevation;
    pImpl->mMaximumElevation = elevation;
    pImpl->mHaveTopography = true;
}

/// Have toppography?
bool Constant::haveTopography() const noexcept
{
   return pImpl->mHaveTopography;
}

/// Evaluate
double Constant::evaluate(const double, const double) const
{
    return evaluate(0, 0, nullptr, nullptr);
} 

/// Evaluate with derivatives
double Constant::evaluate(const double, const double,
                          double *dElevationDx, double *dElevationDy) const
{
    if (!haveTopography()){throw std::runtime_error("Topography not set");}
    if (dElevationDx != nullptr){*dElevationDx = 0;}
    if (dElevationDy != nullptr){*dElevationDy = 0;}
    return pImpl->mElevation;
}

/// Min/max elevation
std::pair<double, double> Constant::getMinimumAndMaximumElevation() const
{
    if (!haveTopography()){throw std::runtime_error("Topography not set");}
    return std::pair {pImpl->mMinimumElevation, pImpl->mMaximumElevation};
}
