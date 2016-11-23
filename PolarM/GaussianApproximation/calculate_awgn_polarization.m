function [ channels ] = calculate_awgn_polarization( llr_vec, n )
% This function calculates polarization

N  = 2^n;
% global bit_reversed; % = bitrevorder((1:N)');

if length(llr_vec) == 1
    channels = llr_vec*ones(1, N);
elseif length(llr_vec) == N
    channels = llr_vec;
else
    disp('Unsupoorted llr, n combinationation');
end

for i=1:n
    c1 = channels(1:2:N);
    c2 = channels(2:2:N);
    channels = [phi_x_inv(1 - (1 - phi_x_table(c1)).*(1 - phi_x_table(c2)))', c1 + c2];
end


end

