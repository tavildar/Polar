function [ mean_llr ] = get_bpsk_llr_for_capacity( const_capacity )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


mean_llr = 0 * const_capacity;

load('bpsk_cap.mat');

for i = 1:size(const_capacity,1)
    for j = 1:size(const_capacity,2)
        for k = 1:size(const_capacity,3)
            for snr_index = 1 : length(snr_vec_db)
                if capacity(snr_index) >= const_capacity(i, j, k)
                    break;
                end
            end
            mean_llr(i, j, k) = 4 * 10^(snr_vec_db(snr_index)/10);
        end
    end
end
end

