#ifndef ULOCATOR_H5IO_HPP
#define ULOCATOR_H5IO_HPP
#include <string>
#include <vector>
#include <hdf5.h>
namespace
{
herr_t h5read(hid_t dataSet, hid_t memSpace, hid_t dataSpace,
              std::vector<float> *result)
{
    return H5Dread(dataSet, H5T_NATIVE_FLOAT, memSpace, dataSpace,
                   H5P_DEFAULT, result->data());
}
herr_t h5read(hid_t dataSet, hid_t memSpace, hid_t dataSpace,
              std::vector<double> *result)
{
    return H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace, dataSpace,
                   H5P_DEFAULT, result->data());
}
template<typename T>
void readDataset(const hid_t mGroup,
                 const std::string &dataSetName,
                 std::vector<T> *result,
                 std::vector<size_t> *dimsOut = nullptr)
{
    auto dataSet = H5Dopen2(mGroup, dataSetName.c_str(), H5P_DEFAULT);
    auto dataSpace = H5Dget_space(dataSet);
    auto rank = H5Sget_simple_extent_ndims(dataSpace);
    std::vector<hsize_t> dims(rank);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    hsize_t length = 1;
    if (dimsOut != nullptr)
    {
        dimsOut->clear();
        dimsOut->reserve(dims.size());
    }
    for (const auto &di : dims)
    {
        length = length*di;
        if (dimsOut != nullptr){dimsOut->push_back(static_cast<size_t> (di));}
    }
    try
    {
        result->resize(length, 0);
    }
    catch (...)
    {
        H5Sclose(dataSpace);
        H5Dclose(dataSet);
        throw std::runtime_error("Failed to resize result");
    }
    // Read the data
    auto memSpace = H5Screate_simple(rank, dims.data(), nullptr);
    auto status = ::h5read(dataSet, memSpace, dataSpace, result);
    H5Sclose(memSpace);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    if (status != 0)
    {
        result->clear();
        throw std::runtime_error("Failed to read dataset: " + dataSetName);
    }
}
}
#endif
