#include "shapelet_transform.hpp"

TimeSeries::TimeSeries(const std::vector<float>& values, const int class_id)
        : values(values), class_id(class_id)
{
    //empty
}


Shapelet::Shapelet(TimeSeries& ts, int start_index, int length)
:   ts(&ts), class_id(ts.get_class_id()), start_index(start_index), length(length)
{
    //empty
}


// compare if two shapelets are self similar
bool Shapelet::operator==(Shapelet& s){
    // compare if both shapelets are from the same time series and
    // if they share any time series elements
    return (ts == s.get_TimeSeries()) and
        ((this->get_start_position() >= s.get_start_position() and this->get_start_position() < s.get_final_posiiton()) or
        (s.get_start_position() >= this->get_start_position() and s.get_start_position() < this->get_final_posiiton()));
}


ShapeletSelection::ShapeletSelection(std::vector<TimeSeries> ts_list, const unsigned int k_best,
                                    const int min, const int max)
: ts_list(ts_list), k(k_best), min(min), max(max)
{
    if(min > max or min < 0){
        throw std::invalid_argument("Error initializing ShapeletSelection:  \
        Min greater than Max!");
    }
    if(ts_list.size() <= 2){
        throw std::invalid_argument("Error initializing ShapeletSelection:  \
        Number of time series must be greater than 2!");
    }
}


// initialize reading from a dataset file
// headless should be false if no csv file is present
ShapeletSelection::ShapeletSelection(std::string input_filepath, const int k_best,
                                    const int min, const int max)
: ts_list(read_dataset(input_filepath)), k(k_best), min(min), max(max)
{
    if(min > max or min < 0){
        std::cerr << "Min greater than max!" << std::endl;
        exit(-1);
    }
    if(ts_list.size() <= 2){
        std::cerr << "Number of time series must be greater than 2" << std::endl;
        exit(-2);
    }
}


void ShapeletSelection::select_best()
{
    int num_shapelets=0; //defines the number of shapelets of size l
    // holds the value of each distance calculations.
    // Cleared after calculation is complete and reused;
    std::vector<float> shapelet_distances;
    std::vector<Shapelet> candidate_list;

    // for every time series in ts_list
    for(TimeSeries& current_ts: ts_list)
    {
        // for each length between min and max
        for(int l=min; l <= max; l++){
            num_shapelets = current_ts.size() - l + 1;  // number of shapelets of size l
            //for each position of a given length
            for(int position=0; position < num_shapelets; position++){
                Shapelet shapelet_candidate(current_ts, position, l);

                // Calculate distances from current shapelet candidate to each time series in T,
                for(TimeSeries& compared_ts : ts_list){
                   shapelet_distances.push_back(shapelet_ts_distance(shapelet_candidate, compared_ts));
                }

                shapelet_candidate.set_quality(bin_f_statistic(shapelet_distances));
                candidate_list.push_back(shapelet_candidate);

                //clear distances vector
                shapelet_distances.clear();
            }
        } // Here all shapelets and qualities from current_ts should be calculated
        // sort shapelets in descending order by quality
        sort(candidate_list.begin(), candidate_list.end(), std::greater<Shapelet>());
        remove_self_similars(candidate_list);
        // merge found shapelets with current best
        best_shapelets.insert(best_shapelets.end(), candidate_list.begin(), candidate_list.end());
        sort(best_shapelets.begin(), best_shapelets.end(), std::greater<Shapelet>());
        // trim the list to only keep the K best
        if(best_shapelets.size() > k){
            best_shapelets.erase(best_shapelets.begin()+k, best_shapelets.end());
        }

        candidate_list.clear();
    }
}


// calculates eucledian distance between pivot and target.
// Modifies minimun distance if a new smaller value is found, else
// does not modify anything
void ShapeletSelection::euclidean_distance(const std::vector<float>& pivot,
                                          const std::vector<float>& target,
                                          float& curr_min)
{
   float total_dist = 0;

    for(int i=0; i < pivot.size(); i++){
        total_dist += pow(pivot.at(i) - target.at(i), 2);
        // early abandon
        if(total_dist >= curr_min) return;
    }
    // new value is smaller than current, update current minimun;
    curr_min = total_dist;
}


// normalizes float vector std::vector<float> v (overwrites the vector)
void ShapeletSelection::zscore_normalization(std::vector<float>& vec)
{
    float mean=0, std, diff_sum=0;
    for(float v : vec){
        mean+=v;
    }
    mean /= (float) vec.size();

    for(float v : vec){
        diff_sum+= pow(v - mean, 2);
    }
    diff_sum /= (float) (vec.size() - 1);

    std = sqrt(diff_sum);

    // special case, std results in zero, normalizing everything to zero
    if(std == 0){
        fill(vec.begin(), vec.end(), 0.0);
    }
    else{
    // normalize vector
        for(float& v : vec){
            v = (v - mean)/std;
        }
    }
}


// calculates the distance between a single ts to an entire timeseries
// returns the smallest distance between this shapelet and all same size shapelets
// in the time series
float ShapeletSelection::shapelet_ts_distance(Shapelet& pivot, TimeSeries& curr_ts)
{
    float min_distance = INFINITY;
    std::vector<float> pivot_norm, target_norm;
    const int num_shapelets = curr_ts.size() - pivot.get_length() + 1;

    // normalize pivot only ONCE
    pivot_norm = pivot.get_values_vector();
    zscore_normalization(pivot_norm);

    for(int i=0; i< num_shapelets; i++){
        // copy only the current shapelet values necessary
        target_norm = std::vector<float>(curr_ts.get_values().begin() + i,
                                    curr_ts.get_values().begin() + i + pivot.get_length());
        zscore_normalization(target_norm);  // normalize target shaplet

        euclidean_distance(pivot_norm, target_norm, min_distance);
    }
    return min_distance;
}


// this code is awful and i feel bad
float ShapeletSelection::bin_f_statistic(std::vector<float>& distances)
{
    float class_zero_sum=0, class_one_sum=0, total_dist_sum;
    float class_zero_avg, class_one_avg, total_dist_avg;
    float numerator_sum, denominator_sum=0;
    int class_zero_ts_num=0, class_one_ts_num=0;
    // Count the number of time series in each class and compute the sum of distaces for each class
    for(int i=0; i < ts_list.size(); i++){
        if(ts_list.at(i).get_class_id() == 0){
            class_zero_sum += distances.at(i);
            class_zero_ts_num++;
        }
        else if (ts_list.at(i).get_class_id() == 1){
            class_one_sum += distances.at(i);
            class_one_ts_num++;
        }
        else{
            std::cout << "Class is not binary" << std::endl;
            exit(-1);
        }
    }

    //Claculate the average values for each class and for the entire distance
    class_zero_avg = class_zero_sum/class_zero_ts_num;
    class_one_avg = class_one_sum/class_one_ts_num;
    total_dist_sum = class_zero_sum + class_one_sum;
    total_dist_avg = total_dist_sum/ts_list.size();

    // f-stat formula numerator
    numerator_sum = pow(class_zero_avg - total_dist_avg, 2) + pow(class_one_avg - total_dist_avg, 2);

    // Calculates the sums in f-stat formula denominator
    for(int i=0; i < ts_list.size(); i++)
    {
        if (ts_list.at(i).get_class_id() == 0){
            denominator_sum += pow(distances.at(i) - class_zero_avg, 2);
        }
        else {      // non binary classes already verified
            denominator_sum += pow(distances.at(i) - class_one_avg, 2);
        }
    }
    return numerator_sum/(denominator_sum/(ts_list.size() -2));
}


void ShapeletSelection::remove_self_similars(std::vector<Shapelet>& shapelet_list)
{
    //new list containing non self similar shapelets
    std::vector<Shapelet> temp_list;
    bool self_similar=false;

    //first shapelet is always included, given a sorted list
    temp_list.push_back(shapelet_list.at(0));
    for(int i=1; i < shapelet_list.size(); i++){
        //check every other shapelet in the list from begging up to itself for self similars
        self_similar = false;
        for(int j=0; j < i; j++){
            if(shapelet_list.at(i) == shapelet_list.at(j)){
                self_similar = true;
                break;
            }
        }
        // not self similar shapelets are added to temp_list
        if(!self_similar){
            temp_list.push_back(shapelet_list.at(i));
            // stop checking for self similars if k new best shapelets have been added
            if(temp_list.size() >= k){
                break;
            }
        }
    }

    //copy new list to old one
    shapelet_list = temp_list;
}


// read csv file header must have format
// <number_ts> <number_elements>
// SEPARATED BY SPACE
// each following line is a ts, elements separated by spaces, and
// the last element is 0 or 1 representing class
std::vector<TimeSeries> ShapeletSelection::read_dataset(std::string input_filepath)
{
    unsigned long num_ts, num_elements; // number of timeseries and elements per ts
    std::vector<TimeSeries> ts_list;
    std::string line;
    std::vector<float> values;
    int class_id;
    float val;

    std::ifstream input_file(input_filepath);
    if(input_file.fail()){
        std::cerr << "Could not open file: " << input_filepath << " for reading!" << std::endl;
        exit(-1);
    }

    //reads and parse header
    std::getline(input_file, line);
    std::stringstream _s(line);
    _s >> num_ts;
    _s >> num_elements;

    if(num_ts <= 0 or num_elements <= 0){
        std::cerr << "Error reading csv header!" << std::endl;
        exit(-1);
    }

    // Read data, line by line
    while(std::getline(input_file, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);

        // Extract each integer
        while(ss >> val){
            values.push_back(val);
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
        }
        // last element is class
        class_id = (int) values.back();
        values.pop_back(); // remove class id
        // check values vector has the correct amount of elements
        if(values.size() != num_elements){
            std::cerr << "Error reading timeseries [" << ts_list.size() + 1 <<
            "] with " << values.size() << " elements, expected: " << num_elements << std::endl;
            exit(-1);
        }
        // add timeseries to list
        ts_list.push_back(TimeSeries(values, class_id));
        // clear values vector
        values.clear();
    }
    // check the number of ts read
    if(ts_list.size() != num_ts){
        std::cerr << "Error: TS read: " << ts_list.size() << " TS expected: " <<
        num_ts << std::endl;
    }
    input_file.close();

    return ts_list;
}


void ShapeletSelection::write_best_shapelets(std::string output_filepath) const
{
    std::ofstream output_file(output_filepath);
    if(output_file.fail()){
        std::cerr << "Error opening output file!" << std::endl;
        exit(-1);
    }

    for(auto s : best_shapelets){
        for(auto v : s.get_values_vector()){
            output_file << v << ",";
        }
        output_file << s.get_quality() << "\n";
    }

    output_file.close();
}

//------------------------
// OPEN CL IMPLEMENTATION
//------------------------

// lets user choose a platform based on availble platforms
void ShapeletSelectionCL::interactive_set_platform_device()
{
    std::vector<cl::Platform> platforms;    // platform query

    // get platforms
    cl::Platform::get(&platforms);

    std::cout << "Select the platform:" << std::endl;

    for(int i=0; i < platforms.size(); i++){
        std::cout << "[" << i << "] " << platforms.at(i).getInfo<CL_PLATFORM_NAME>() << std::endl;
    }

    int index = -1;
    while(index < 0 or index >= platforms.size()){
        std::cin >> index;
    }

    platform = platforms.at(index);
    std::cout << "Selected platform: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

    // get devices from platform
    platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);

    std::cout << "Select device: " << std::endl;

    for(int i=0; i < devices.size(); i++){
        std::cout << "[" << i << "] " << devices[i].getInfo<CL_DEVICE_NAME>() << std::endl;
    }

    index = -1;
    while(index < 0 or index >= devices.size()){
        std::cin >> index;
    }

    device = devices.at(index);
    std::cout << "Selected device: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
}


// selects the first platform and device found
void ShapeletSelectionCL::set_default_platform_device()
{
    std::vector<cl::Platform> platforms;    // platform query

    cl::Platform::get(&platforms);
    platform = platforms.at(0);

    platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
    device = devices.at(0);
}


// creates context, command queue, buffers and kernel
void ShapeletSelectionCL::init_context()
{
    context = cl::Context(devices);

    if(enable_profiling){
        queue = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE);
    }
    else {
        queue = cl::CommandQueue(context, device);
    }

    // create and compile kernel
    std::ifstream sourceFile(KERNEL_FILEPATH);
    if(sourceFile.fail()){
        std::cerr << "Could not find kernel file" << std::endl;
        exit(-1);
    }

    std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),
                        (std::istreambuf_iterator < char > ()));

    cl::Program::Sources source(1,
                            std::make_pair(sourceCode.c_str(),
                            sourceCode.length() + 1));

    // Make program from the source code
    cl::Program program(context, source);

    try
    {
        program.build(devices);
    }
    catch (cl::Error& e)
    {
        if (e.err() == CL_BUILD_PROGRAM_FAILURE)
        {
            for (cl::Device dev : devices)
            {
                // Check the build status
                cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev);
                if (status != CL_BUILD_ERROR)  continue;

                // Get the build log
                std::string name     = dev.getInfo<CL_DEVICE_NAME>();
                std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev);
                std::cerr << "Build log for " << name << ":" << std::endl
                            << buildlog << std::endl;
            }
        }
        else
        {
            throw;
        }
    }
    // Get the build log for device
    build_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device); 

    // Make kernel
    kernel = cl::Kernel(program, "shapelet_distance");

    // create new buffers
    // pivot is constant memory access
    // each buffer is created with the maximum memory required predicted
    buf_norm_pivot = cl::Buffer(context, CL_MEM_READ_ONLY, max * sizeof(float));
    // entire time series to access target shapelets. Each TS must have the same size
    // so it is ok to set this buffer size to be the size of the first time series
    buf_target_shapelet = cl::Buffer(context, CL_MEM_READ_ONLY, ts_list[0].size() * sizeof(float));
    // heap memory area to store normalized pivot shapelets
    // the maximum number of shapelets is ts.size() - min + 1, and its maximum length is max
    buf_norm_target = cl::Buffer(context, CL_MEM_READ_WRITE,  (ts_list[0].size() - min + 1) * max * sizeof(float));
    // stores distance results
    buf_distances = cl::Buffer(context, CL_MEM_WRITE_ONLY, (ts_list[0].size() - min + 1) * sizeof(float));

    // set buffers as kernel args
    // the shapelet length must be changed each kernel execution
    kernel.setArg(0, buf_norm_pivot);
    kernel.setArg(1, buf_target_shapelet);
    kernel.setArg(2, buf_norm_target);
    kernel.setArg(3, buf_distances);

}


ShapeletSelectionCL::ShapeletSelectionCL(const std::string input_filepath, const int k_best,
    const int min, const int max, const bool interactive_selection, const bool enable_profiling)
: ShapeletSelection(input_filepath, k_best, min, max), enable_profiling(enable_profiling)
{
    try{
        if(interactive_selection){
            interactive_set_platform_device();
        }
        else{
            set_default_platform_device();
        }

        init_context();
    }
    catch(cl::Error& error){
        std::cout << error.what() << "(" << error.err() << ")" << std::endl;
        exit(-1);
    }
}


float ShapeletSelectionCL::shapelet_ts_distance(Shapelet& pivot, TimeSeries& curr_ts)
{
    std::vector<float> pivot_norm;
    const int num_shapelets = curr_ts.size() - pivot.get_length() + 1;
    float *distances = new float[num_shapelets]; // result distances
    float min_distance;

    // normalize pivot shapelet only once on host
    pivot_norm = pivot.get_values_vector();
    zscore_normalization(pivot_norm);

    // launch kernels from here
    // each work-item will calculate the distance between pivot and a target shapelet
    try
    {
        // write the new data to buffer arguments
        queue.enqueueWriteBuffer(buf_norm_pivot, CL_TRUE, 0,  pivot_norm.size() * sizeof(float), pivot_norm.data());
        queue.enqueueWriteBuffer(buf_target_shapelet, CL_TRUE, 0,  curr_ts.size() * sizeof(float), curr_ts.get_values().data());

        // set shapelet length
        kernel.setArg(4, (int) pivot_norm.size());

        // the number of work items is equal to the number of shapelets to compare
        cl::NDRange global(num_shapelets);

        queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, 
                    cl::NullRange, NULL, &event);

        if(enable_profiling){
            event.wait();
            // time necessary to issue the kernel
            total_dispatch_time += event.getProfilingInfo<CL_PROFILING_COMMAND_SUBMIT>() -
                                    event.getProfilingInfo<CL_PROFILING_COMMAND_QUEUED>();

            // time executing the kernel
            total_exec_time += event.getProfilingInfo<CL_PROFILING_COMMAND_END>() -
                            event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
        }
        queue.finish();

        // write values to a distance array
        queue.enqueueReadBuffer(buf_distances, CL_TRUE, 0, num_shapelets * sizeof(float), distances);
    }
    catch(cl::Error& error)
    {
        std::cerr << "Error during kernel execution! " << error.what() <<
        " [" << error.err() << "] " << std::endl;
        exit(-1);
    }

    // find the min distance in array
    min_distance = *std::min_element(distances, distances + num_shapelets);

    delete[] distances;

    return min_distance;
}


void ShapeletSelectionCL::print_device_info() const
{
    std::cout << "Platform: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl
    << "Device: " << device.getInfo<CL_DEVICE_NAME>() << std::endl
    << "Type: " << device.getInfo<CL_DEVICE_TYPE>() << std::endl;
}


void ShapeletSelectionCL::print_build_log() const
{
    std::cout << "Device: " << device.getInfo<CL_DEVICE_NAME>() << std::endl
    << build_log << std::endl;
}
