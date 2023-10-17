#ifndef ULOCATOR_QUARRY_HPP
#define ULOCATOR_QUARRY_HPP
#include <memory>
namespace ULocator::Position
{
 class WGS84;
}
namespace ULocator
{
class Quarry
{
public:
    /// @brief Constructor.
    Quarry();
    /// @brief Copy constructor.
    Quarry(const Quarry &quarry);
    /// @brief Move constructor.
    Quarry(Quarry &&quarry) noexcept;
    /// @brief Creates a quarry with a given location and name.
    Quarry(const Position::WGS84 &position, const std::string &name);

    void setGeographicPosition(const Position::WGS84 &position);
    [[nodiscard]] Position::WGS84 getGeographicPosition() const;
    [[nodiscard]] bool haveGeographicPosition() const noexcept;

    /// @brief Sets the quarry name.
    void setName(const std::string &name);
    /// @result The quarry name.
    [[nodiscard]] std::string getName() const noexcept;
    
    /// @brief Resets the class. 
    void clear() noexcept;
    /// @brief Destructor
    ~Quarry();

    Quarry& operator=(const Quarry &quarry);
    Quarry& operator=(Quarry &&quarry) noexcept;
private:
    class QuarryImpl;
    std::unique_ptr<QuarryImpl> pImpl;
};
}
#endif
