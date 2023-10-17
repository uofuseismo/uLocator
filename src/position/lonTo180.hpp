#ifndef PRIVATE_TO180_HPP
#define PRIVATE_TO180_HPP
#include <cmath>
#include <limits>
#ifndef NDEBUG
#include <cassert>
#endif
namespace
{
template<typename T>
[[nodiscard]]
T lonTo180(const T lonIn)
{
    auto lon = lonIn;
    if (lon < -180)
    {
        for (int k = 0; k < std::numeric_limits<int>::max(); ++k)
        {
            auto tempLon = lon + k*360;
            if (tempLon >= -180)
            {   
                lon = tempLon;
                break;
            }
        }
    }   
    if (lon >= 180)
    {
        for (int k = 0; k < std::numeric_limits<int>::max(); ++k)
        {
            auto tempLon = lon - k*360;
            if (tempLon < 180)
            {   
                lon = tempLon;
                break;
            }
        }
    }   
#ifndef NDEBUG 
    assert(lon >= -180 && lon < 180);
#endif
    return lon;
}
}
#endif
