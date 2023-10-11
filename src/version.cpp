#include <string>
#include "uLocator/version.hpp"

using namespace ULocator;

int Version::getMajor() noexcept
{
    return ULOCATOR_MAJOR;
}

int Version::getMinor() noexcept
{
    return ULOCATOR_MINOR;
}

int Version::getPatch() noexcept
{
    return ULOCATOR_PATCH;
}

bool Version::isAtLeast(const int major, const int minor,
                        const int patch) noexcept
{
    if (ULOCATOR_MAJOR < major){return false;}
    if (ULOCATOR_MAJOR > major){return true;}
    if (ULOCATOR_MINOR < minor){return false;}
    if (ULOCATOR_MINOR > minor){return true;}
    if (ULOCATOR_PATCH < patch){return false;}
    return true;
}

std::string Version::getVersion() noexcept
{
    std::string version(ULOCATOR_VERSION);
    return version;
}
