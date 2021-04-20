#pragma once

#define __CL_ENABLE_EXCEPTIONS
#define CL_TARGET_OPENCL_VERSION 120

#include <CL/cl.hpp>
#include <algorithm> // for std::sort
#include <cmath>     // for pow()
#include <fstream>
#include <iostream>
#include <sstream> // for csv file reading
#include <string>
#include <vector>

// TODO: move member functions from hpp to cpp
// time series is a vector with length and associated class
class TimeSeries {
private:
    const std::vector<float> values;
    const int class_id; // stores class identifier
public:
    TimeSeries(const std::vector<float> &values, const int class_id);

    // getters
    const std::vector<float> &get_values() const { return values; }

    int get_class_id() const { return class_id; }

    size_t size() const { return values.size(); }

    float operator[](size_t i) const { return values[i]; }
};

// TODO: move member functions from hpp to cpp
// IMPORTANT: a shapelet refers to a TimeSeries for its values
// Thus, the given TS must not be destroyed afterward
class Shapelet {
private:
    TimeSeries *ts; // used only to infer if two shapelets are self similar and get shapelet values
    int class_id;
    int start_index, length;
    float quality = 0; // not every shapelet has quality, and is only set after being calculated

public:
    Shapelet(TimeSeries &ts, int start_index, int length);

    // gettes and setters
    void set_quality(float q) { quality = q; }

    int get_class_id() const { return class_id; }

    float get_quality() const { return quality; }

    // returns a copy of shapelet's values
    std::vector<float> get_values_vector() const {
        return std::vector<float>(ts->get_values().begin() + start_index,
                                  ts->get_values().begin() + get_final_posiiton());
    }

    int get_length() { return length; }

    // this is dumb but necessary for identifying the original time series. TODO: find better
    // solution
    const TimeSeries *get_TimeSeries() const { return ts; }

    int get_start_position() const {
        return start_index;
    } // start posiiton in relation to time series

    int get_final_posiiton() const {
        return start_index + length;
    } // end position in relation to time series

    // operator used to sort by quality
    bool operator>(const Shapelet &s) const { return quality > s.get_quality(); }

    // checks if two shapelets are self similar
    // no need to be implemented as operator, I just wanted to use it this way
    bool operator==(const Shapelet &s);
};

class ShapeletSelection {
protected:
    std::vector<TimeSeries> ts_list;      // list cointaining every time series
    std::vector<Shapelet> best_shapelets; // list with the K best shapelets
    const unsigned int k;                 // number of best shapelets to save
    int min, max;                         // min and max values for each shapelet

    // internal methods
    void euclidean_distance(const std::vector<float> &pivot, const std::vector<float> &target,
                            float &curr_min);
    void zscore_normalization(std::vector<float> &v);
    float bin_f_statistic(std::vector<float> &distances);
    void remove_self_similars(std::vector<Shapelet> &shapelet_list);

    // overriden by class ShapeletSelectionCL
    virtual float shapelet_ts_distance(std::vector<float> &pivot, TimeSeries &curr_ts);

public:
    ShapeletSelection(const std::vector<TimeSeries> &ts_list, const unsigned int k_best,
                      const int min, const int max);

    // initialize reading from a csv file
    ShapeletSelection(const std::string &input_filepath, const int k_best, const int min,
                      const int max);

    // virtual destructor
    virtual ~ShapeletSelection() = default;

    // runs shapelet selection to find the k best shapelets
    void select_best();

    // getter and setters
    const std::vector<Shapelet> &get_best_shapelets() const { return best_shapelets; }

    // read csv file
    // HEADER CONTAINTS: <number of ts> <ts elements>
    // header is separated by space
    // each line is a time series separated by comma, last value is 0 or 1 representing class
    std::vector<TimeSeries> read_dataset(const std::string &input_filepath);

    // writes the best shapelets to file
    // comma separated, no header
    void write_best_shapelets(const std::string &output_filepath) const;
};

// OpenCL accelerated implementation of ShapeletSelection algoritm
class ShapeletSelectionCL : public ShapeletSelection {
private:
    cl::Platform platform;
    cl::Device device;
    cl::Context context;
    cl::CommandQueue queue;
    cl::Event event;                 // event exclusive for profiling
    std::vector<cl::Device> devices; // devices query
    cl::Kernel kernel;
    cl::Buffer buf_norm_pivot, buf_target_shapelet, buf_norm_target, buf_distances;
    std::string build_log; // kernel build log
    // profiling variables
    const bool enable_profiling;      // true to profile kernel execution time
    cl_ulong total_exec_time = 0;     // total execution time for kernels
    cl_ulong total_dispatch_time = 0; // total time dispatching the kernel

    const char *KERNEL_FILEPATH = "kernels/shapelet_distance_kernel.cl";
    const char *KERNEL_FUNCTION_NAME = "shapelet_distance";

    // user interactive select for platform and device
    void interactive_set_platform_device();

    // selects the first platform an device found
    void set_default_platform_device();

    // init context, command queue, buffers and kernel
    // perf_prof: adds kernel performance profiling to kernel execution
    // prints information to stdout
    void init_context();

protected:
    // shapelet_ts_distance is overriden to a accelerated implementation
    float shapelet_ts_distance(std::vector<float> &pivot_norm, TimeSeries &curr_ts) override;

public:
    // interactive_selection: selects the first platform and device found from query
    // otherwise user input must select platform an device
    // enable_profiling: adds profiling to CL enqueue and
    // prints kernel perfomance information to stdout
    ShapeletSelectionCL(std::string input_filepath, const int k_best, const int min, const int max,
                        const bool interactive_selection = true,
                        const bool enable_profiling = true);

    // print platform and device information
    void print_device_info() const;

    // TODO: IMPLEMENT
    // prints compiler kernel build info for current selected device
    void print_build_log() const;

    cl_ulong get_total_exec_time() const { return total_exec_time; }
    cl_ulong get_total_dispatch_time() const { return total_exec_time; }
};

