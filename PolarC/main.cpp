#include <iostream>
#include <iomanip>      // std::setprecision
#include <vector>
#include "PolarCode.h"

int main(int argc, char* argv[])
{

    uint8_t n = 11;
    uint16_t info_length = (uint16_t) (1 << (n - 1));
    uint16_t crc_size = 0;

    double design_epsilon = 0.32;

    PolarCode polar_code(n, info_length, design_epsilon, crc_size);

    double ebno_log_min = 1.00;
    double ebno_log_max = 2.01;
    double ebno_log_increment = 0.25;
    std::vector<double> ebno_vec;

    for (double ebno_log = ebno_log_min; ebno_log <= ebno_log_max; ebno_log += ebno_log_increment)
        ebno_vec.push_back(ebno_log);

    std::vector<uint8_t> list_size_vec(0);
    list_size_vec.push_back(1);
    list_size_vec.push_back(2);
    list_size_vec.push_back(4);
    list_size_vec.push_back(8);
    list_size_vec.push_back(32);

    std::vector<std::vector<double>> bler = polar_code.get_bler_quick(ebno_vec, list_size_vec);

    for (unsigned i_ebno = 0; i_ebno < ebno_vec.size(); ++i_ebno) {
        std::cout << std::fixed  << std::setprecision(3) << ebno_vec.at(i_ebno) << "\t \t";
        for (unsigned i_list = 0; i_list < list_size_vec.size(); ++i_list) {
            std::cout << std::fixed  << std::setprecision(6) << bler.at(i_list).at(i_ebno) << "\t";
        }
        std::cout << std::endl;
    }

    return 0;
}
