#ifndef ULOCATOR_POSITION_KNOWN_LOCAL_LOCATION_HPP
#define ULOCATOR_POSITION_KNOWN_LOCAL_LOCATION_HPP
#include <memory>
namespace ULocator::Position
{
/// @class IKnownLocalLocation "knownLocalLocation.hpp" "uLocator/position/knownLocalLocation.hpp"
/// @brief A known location is simply an (x,y,z) location in the local coordinates.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class IKnownLocalLocation
{
public:
    /// @result The local x coordinates in meters.
    [[nodiscard]] virtual double x() const = 0;
    /// @result The local y coordinates in meters.
    [[nodiscard]] virtual double y() const = 0;
    /// @result The local z coordinates in meters.  This increases + down.
    [[nodiscard]] virtual double z() const = 0;
    /// @result A copy of the base class.
    [[nodiscard]] virtual std::unique_ptr<IKnownLocalLocation> clone() const = 0;
    /// @brief Destructor
    virtual ~IKnownLocalLocation() = default;
};
}
#endif
