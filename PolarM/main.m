
%main m
clear all;
n = 11;
block_length = 2^n;
info_length = 2^(n-2);
crc_size = 0;
epsilon = 0.32;
list_vec = [1, 4];
ebno_vec = 1:0.25:2;
polar_code = PolarCode(block_length, info_length, epsilon, crc_size);
[bler, ~] = polar_code.get_bler_quick(ebno_vec, list_vec);
[ebno_vec', bler]
