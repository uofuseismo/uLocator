#include <uLocator/version.hpp>
#include <pybind11/pybind11.h>
#include "dataStructures.hpp"
#include "position.hpp"
#include "topography.hpp"
#include "rayTracers.hpp"
//#include "broadcasts.hpp"
//#include "services.hpp"

PYBIND11_MODULE(pyulocator, m)
{
    m.attr("__version__") = ULOCATOR_VERSION;
    m.attr("__name__") = "pyulocator";
    m.attr("__doc__") = "A Python interface to the Univeristy of Utah Seismograph Stations Location library.";

    ULocator::Python::initialize(m);
    ULocator::Python::Position::initialize(m);
    ULocator::Python::Topography::initialize(m);
    ULocator::Python::RayTracers::initialize(m);
}

