function [ c ] = get_bpsk_cap( snr_db )
% assumes e[x^2] = 1;

x_vec = [-1, 1];
n_0 = 1/2*10^(-snr_db/10);
delta_y = sqrt(n_0) * 0.001;
max_value = min(10000, max(x_vec) + 3 + 3*sqrt(n_0));
y_vec = -max_value: delta_y : max_value;
p_y = zeros(length(y_vec), 1);
for y_index = 1 : length(y_vec)
    y = y_vec(y_index);
    for x_index = 1 : length(x_vec) 
        x = x_vec(x_index);
        p_y(y_index) = p_y(y_index) + exp(-(y - x).^2/2/n_0)/sqrt(2*pi*n_0) * 0.5;
    end
end

p_y = p_y / (sum(p_y)*delta_y);

h = 0;
for y_index = 1 : length(y_vec)
    if p_y(y_index) > 0
        h = h + (-log2(p_y(y_index))) * p_y(y_index) * delta_y;
    end
end
c = h - 0.5 * (1 + log(2*pi*n_0))/log(2) ;


