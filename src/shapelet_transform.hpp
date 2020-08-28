#ifndef _SHAPELET_TRANSFORM_H
#define _SHAPELET_TRANSFORM_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>   // for csv file reading
#include <algorithm> //for std::sort
#include <cmath>    // for pow()

// TODO:
// Add exceptions for errors


// TODO: move member functions from hpp to cpp
// time series is a vector with length and associated class
class TimeSeries
{
    private:
        const std::vector<float> values;
        const int class_id;    // stores class identifier
    public:
        TimeSeries(const std::vector<float> values, const int class_id)
            : values(values), class_id(class_id)
        {
            //empty constructor body
        }

        // getters
        const std::vector<float>& get_values() const { return values; }

        int get_class_id() const { return class_id; }

        int size() const{ return values.size(); }

        const float operator[](int i) const { return values[i]; }
};


// TODO: move member functions from hpp to cpp
// IMPORTANT: a shapelet refers to a TimeSeries for its values
// Thus, the given TS must not be destroyed afterward
class Shapelet
{
    private:    
        TimeSeries* ts;     // used only to infer if two shapelets are self similar and get shapelet values
        int class_id;
        int start_index, length;
        float quality = 0;    //not every shapelet has quality, and is only set after being calculated
    
    public:
        Shapelet(TimeSeries& ts, int start_index, int length)
        :   ts(&ts), class_id(ts.get_class_id()), start_index(start_index), length(length)
        {
            //empty
        }

        // gettes and setters
        void set_quality(float q) { quality = q; }

        int get_class_id() const { return class_id ;}

        float get_quality() const { return quality; }

        // returns a copy of shapelet's values
        std::vector<float> get_values_vector() const 
        { 
            return std::vector<float>(ts->get_values().begin()+start_index, 
                                        ts->get_values().begin()+get_final_posiiton()); 
        }

        int get_length(){ return length; }

        // this is dumb but necessary for identifying the original time series. TODO: find better solution
        const TimeSeries* get_TimeSeries() const { return ts; }

        int get_start_position() const { return start_index; }  // start posiiton in relation to time series
        
        int get_final_posiiton() const {return start_index + length; } // end position in relation to time series
        
        // operator used to sort by quality
        bool operator>(const Shapelet& s) const { return quality > s.get_quality(); }

        // checks if two shapelets are self similar
        bool operator==(Shapelet& s){
            // compare if both references are the same address
            if(ts != s.get_TimeSeries()){
                return false;
            }
            if( (this->get_start_position() >= s.get_start_position() and this->get_start_position() < s.get_final_posiiton()) or
                (s.get_start_position() >= this->get_start_position() and s.get_start_position() < this->get_final_posiiton()))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
};


class ShapeletSelection
{
    private:
        std::vector<TimeSeries> ts_list; //list cointaining every time series
        const unsigned int k;                // number of best shapelets to save
        int min, max;         // min and max values for each shapelet
        std::vector<Shapelet> best_shapelets; // list with the K best shapelets

        //internal methods
        void euclidean_distance(std::vector<float>& pivot, std::vector<float>& target, float& curr_min);
        void zscore_normalization(std::vector<float>& v);
        float shapelet_ts_distance(Shapelet& pivot, TimeSeries& curr_ts);
        float bin_f_statistic(std::vector<float>& distances);
        void remove_self_similars(std::vector<Shapelet>& shapelet_list);
        void merge_shapelets();

    public:
        ShapeletSelection(std::vector<TimeSeries> ts_list, const unsigned int k_best, const int min, const int max);

        // initialize reading from a dataset file
        // headless should be true is csv file is not headless
        ShapeletSelection(std::string input_filepath, const int k_best, const int min, const int max);
        

        // runs shapelet selection to find the k best shapelets
        void select_best();

        // getter and setters
        const std::vector<Shapelet>& get_best_shapelets() const { return best_shapelets; }

        // read csv file 
        // HEADER CONTAINTS: <number of ts> <ts elements>
        // header is separated by space
        // each line is a time series separated by comma, last value is 0 or 1 representing class
        std::vector<TimeSeries> read_dataset(std::string input_filepath); 
        //writes the best shapelets to file
        //comma separated, no header

        void write_best_shapelets(std::string output_filepath) const;
};



#endif