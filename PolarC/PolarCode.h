//
// Created by Saurabh Tavildar on 5/17/16.
//

#ifndef POLARC_POLARCODE_H
#define POLARC_POLARCODE_H


#include <cstdint>
#include <vector>
#include <math.h>
#include <stack>          // std::stack

class PolarCode {


public:

    PolarCode(uint8_t num_layers, uint16_t info_length, double epsilon, uint16_t crc_size) :
            _n(num_layers), _info_length(info_length), _design_epsilon(epsilon),
            _crc_size(crc_size), _llr_based_computation(true)
    {
        _block_length = (uint16_t) (1 << _n);
        _frozen_bits.resize(_block_length);
        _bit_rev_order.resize(_block_length);
        create_bit_rev_order();
        initialize_frozen_bits();
    }

    std::vector<uint8_t> encode(std::vector<uint8_t> info_bits);
    std::vector<uint8_t> decode_scl_p1(std::vector<double> p1, std::vector<double> p0, uint16_t list_size);
    std::vector<uint8_t> decode_scl_llr(std::vector<double> llr, uint16_t list_size);

    std::vector<std::vector<double>> get_bler_quick(std::vector<double> ebno_vec, std::vector<uint8_t> list_size);

private:

    uint8_t _n;
    uint16_t _info_length;
    uint16_t _block_length;
    uint16_t _crc_size;

    double _design_epsilon;

    std::vector<uint8_t> _frozen_bits;
    std::vector<uint16_t> _channel_order_descending;
    std::vector<std::vector<uint8_t>> _crc_matrix;
    std::vector<uint16_t> _bit_rev_order;

    void initialize_frozen_bits();
    void create_bit_rev_order();

    std::vector<uint8_t> decode_scl();
    bool _llr_based_computation;

    std::vector<std::vector<double *>> _arrayPointer_LLR;
    std::vector<double> _pathMetric_LLR;

    uint16_t _list_size;

    std::stack<uint16_t> _inactivePathIndices;
    std::vector<uint16_t > _activePath;
    std::vector<std::vector<double *>> _arrayPointer_P;
    std::vector<std::vector<uint8_t *>> _arrayPointer_C;
    std::vector<uint8_t *> _arrayPointer_Info;
    std::vector<std::vector<uint16_t>> _pathIndexToArrayIndex;
    std::vector<std::stack<uint16_t>> _inactiveArrayIndices;
    std::vector<std::vector<uint16_t>> _arrayReferenceCount;

    void initializeDataStructures();
    uint16_t assignInitialPath();
    uint16_t clonePath(uint16_t);
    void killPath(uint16_t l);

    double * getArrayPointer_P(uint16_t lambda, uint16_t  l);
    double * getArrayPointer_LLR(uint16_t lambda, uint16_t  l);
    uint8_t * getArrayPointer_C(uint16_t lambda, uint16_t  l);

    void recursivelyCalcP(uint16_t lambda, uint16_t phi);
    void recursivelyCalcLLR(uint16_t lambda, uint16_t phi);
    void recursivelyUpdateC(uint16_t lambda, uint16_t phi);

    void continuePaths_FrozenBit(uint16_t phi);
    void continuePaths_UnfrozenBit(uint16_t phi);

    uint16_t findMostProbablePath(bool check_crc);

    bool crc_check(uint8_t * info_bits_padded);

};


#endif //POLARC_POLARCODE_H
