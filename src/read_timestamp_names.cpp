#include <iostream>
#include <vector>
#include "rapidcsv.h"

int main()
{
// "doc" contains all measurements in frame K
rapidcsv::Document tstp_csv("/home/diego/RFS_SLAM/rfsslam/data/timestamps/main_timestamp_list.csv");
std::vector<std::string> list_timestamps_obj;

// rows contains the number of rows (measurements) taken in frame K
unsigned int idx_tstp, rows_tstp_csv = tstp_csv.GetRowCount();
// Iterates over each measurement (x, y, z) and stores them in measurements_k
for(int i = 0; i < rows_tstp_csv ; i++){
    // std::cout << i << std::endl;
    std::string col = tstp_csv.GetRow<std::string>(i)[0];
    list_timestamps_obj.push_back(col);
}

// Print elements of vector col
for(int i = 0; i < rows_tstp_csv ; i++){
    std::cout << list_timestamps_obj[i] << std::endl;
}

return 0;
}