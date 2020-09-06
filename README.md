# ShapeletSelection-OpenCL

 Shapelet Selection Algorithm in C++ based on the [C version](https://github.com/vctrop/shapelet_distance_hardware_accelerator) along with an OpenCL acellerated implementation.

# Usage
A single class is used as interface ShapeletSelection (or ShapeletSelectionCL).
CSV files containing time series can are read, and they MUST follow the format:

number_of_TimeSeries elements_per_TimeSeries

e1, e2, ..., en, class

As such, a header containts the number of TimeSeries contained in the file separeted by a space from the number of elements per TimeSeries. Each TimeSeries must have the exact same number of elements as define by the header.
The following lines are comma separated values, and the last element in each line is its associated class. 

WARNING: currently only 2 classes are accepted for datasets due to the current F statistic implementation. 

A number of DataSets are provided, and an example program extract_shapelets.

# Build:

Just $make . 

An OpenCL library must be installed and available at include and link path. 

This code is written to work with C++11 and OpenCL version 1.2.

# TODO:

* Implement multi class f statistic
* Smarter accleration kernel
* Fix datasets to be compatible with current format
* Reduce memory allocation by making use of C arrays rather than std::vector