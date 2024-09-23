#ifndef ULOCATOR_POSITION_KNOWN_UTAH_QUARRY_HPP
#define ULOCATOR_POSITION_KNOWN_UTAH_QUARRY_HPP
#include <uLocator/position/knownLocalLocation.hpp>
#include <uLocator/position/utahQuarry.hpp>
#include <memory>
#include <vector>
namespace ULocator::Position
{
/// @class KnownUtahQuarry "knownUtahQuarry.hpp" "uLocator/position/knownUtahQuarry.hpp"
/// @brief A known quarry is the (x,y,z) of a Utah quarry in the Utah region.
/// @copyright Ben Baker (University of Utah) distributedd under the MIT license.
class KnownUtahQuarry : public ULocator::Position::IKnownLocalLocation
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.    
    explicit KnownUtahQuarry(const ULocator::Position::UtahQuarry &quarry);
    /// @brief Copy constructor.
    KnownUtahQuarry(const KnownUtahQuarry &quarry);
    /// @brief Move constructor.
    KnownUtahQuarry(KnownUtahQuarry &&quarry) noexcept;
    /// @}

    /// @name Base Class Overrides
    /// @{

    /// @result The local x position of the quarry.
    [[nodiscard]] double x() const override;
    /// @result The local y position of the quarry.
    [[nodiscard]] double y() const override;
    /// @result The local z position of the quarry.
    [[nodiscard]] double z() const override;
    /// @result A copy of this class.
    [[nodiscard]] std::unique_ptr<IKnownLocalLocation> clone() const override;
    /// @}

    /// @name The Unerlying Quarry Information
    /// @{

    /// @result The underlying Utah quarry.
    [[nodiscard]] UtahQuarry getUtahQuarry() const noexcept;
    /// @result A reference to the underlying quarry.
    /// @note This exists for optimization reasons.  You should prefer
    ///       \c getUtahQuarry.
    [[nodiscard]] const UtahQuarry& getUtahQuarryReference() const noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Destructor.
    ~KnownUtahQuarry() override;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] quarry  The quarry to copy to this.
    /// @result A deep copy of the quarry.
    KnownUtahQuarry& operator=(const KnownUtahQuarry &quarry);
    /// @brief Move assignment.
    /// @param[in,out] quarry  The quarry whose memory will be moved to this.
    ///                        On exit, quarry's behavior is undefined.
    /// @result The memory from quarry moved to this.
    KnownUtahQuarry& operator=(KnownUtahQuarry &&quarry) noexcept;
    /// @}

    KnownUtahQuarry() = delete;
private:
    class KnownUtahQuarryImpl;
    std::unique_ptr<KnownUtahQuarryImpl> pImpl;
};
/// @result The known Utah quarries.
std::vector<KnownUtahQuarry> getKnownUtahQuarries();

}
#endif
