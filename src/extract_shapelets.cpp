#include "shapelet_transform.hpp"
#include <iostream>

using namespace std;

// reads a headless csv datafile with a time series followed by a binary class (0 or 1)
// finds the k best shapelets running shapelet selection algorithm
// writes the best shapelets to a csv file
int main(int argc, char *argv[])
{
    int k, min_len, max_len;
    string input_filepath, output_filepath;

    if(argc != 6){
        cerr << "Please use " << argv[0] << 
        " {path_to_dataset} {output_filename} {min_len} {max_len} {k_best}" << endl;
        exit(-1);
    }

    cout << "Shapelet selection algorithm" << endl;
    // parse arguments
    input_filepath = argv[1];
    output_filepath = argv[2];
    min_len = atoi(argv[3]);
    max_len = atoi(argv[4]);
    k = atoi(argv[5]);

    // check if given arguments are valid
    if(k <= 0){
        cerr << "Error: K must be greater than zero!" << endl;
        exit(-1);
    }

    if(min_len < 3){
        cerr << "Error: min_len must be greater than 2!" << endl;
        exit(-1);
    }

    if(min_len > max_len){
        cerr << "min_len must be greater than max_len!" << endl;
        exit(-1);
    }

    // start shapelet selection model
    ShapeletSelectionCL model(input_filepath, k, min_len, max_len);

    cout << "Finding the best " << k << "shapelets from dataset: "
    << input_filepath << " with lengths: " << min_len << " to " << max_len << endl;
    // run shapelet selection algorithm
    model.select_best();
    
    cout << "Shapelet selection complete! Writing best shapelets to: " << output_filepath << endl;
    // write results to file
    model.write_best_shapelets(output_filepath);

    cout << "Sucess!" << endl;
}

