// calculates the distance of a single normalized shapelet
// to an entire ts
// each work-item calculates a distance from a shapelet
// to a shapelet
// results are writen to distances array
__kernel void shapelet_distance(
    __constant float *normalized_pivot_shapelet,    // shapelet must be normalize
    __constant float *target_shapelet,              // contains entire ts. target shapelets are read from this buffer
    __global float *normalized_target_shapelet,     // heap memory to store noramlized shapelets
    __global float *distances,                      // output distance buffer. distance[id] contains output distance
    int shapelet_length
)
{
    float mean=0, diff_sum=0, total_dist=0, std;
    float size = convert_float(shapelet_length);
    int start_position = get_global_id(0);          // start position from target_shapelet
    int offset_position = get_global_id(0) * shapelet_length;   // allocated position in *normalized_target_shapele vector

    // normalize target
    for(int i=0; i < shapelet_length; i++){
        mean += target_shapelet[start_position + i];
    }
    mean /= size;

    for(int i=0; i < shapelet_length; i++){
        diff_sum+= pow(target_shapelet[start_position + i] - mean, 2);
    }
    diff_sum /= (size - 1);

    std = sqrt(diff_sum);

    if(std == 0){
        // std dev is zero, set all values to zeros
        for(int i=0; i < shapelet_length; i++){
            normalized_target_shapelet[offset_position + i] = 0;
        }
    }
    else{
        for(int i=0; i < shapelet_length; i++){
            normalized_target_shapelet[offset_position + i] = (target_shapelet[start_position + i] - mean) / std;
        }
    }

    // calculate distance
    for(int i=0; i < shapelet_length; i++){
        total_dist+=pow(normalized_pivot_shapelet[i] - normalized_target_shapelet[offset_position + i], 2);
    }

    // write distance to output buffer.
    distances[start_position] = total_dist;
}
