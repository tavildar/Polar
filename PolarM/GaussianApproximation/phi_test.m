
%phi 
clear all
close all;
x_vec = 0: 10^(-5): 400;

global minus_log_phi_inv_table;
global min_minus_log_phi;
global max_minus_log_phi;
global increment_minus_log_phi;

phi_vec = (x_vec <  10) .* exp(-0.4527 * x_vec.^(0.86) + 0.0218);
phi_vec = (x_vec >= 10) .* sqrt(pi./(x_vec+0.0001)) .* (1 - 1.4286./(x_vec + 0.0001)) .* exp(-x_vec/4) + phi_vec;
phi_vec = min(phi_vec, 1);
min_minus_log_phi = 0;
max_minus_log_phi = 100;
increment_minus_log_phi = 10^(-5);
minus_log_phi_inv_table = zeros(ceil((max_minus_log_phi - min_minus_log_phi)/increment_minus_log_phi) + 1, 1);

for x_index = 1 : length(x_vec)
   x_value  = x_vec(x_index);
   phi_value = phi_vec(x_index);
   minus_log_phi = -log(phi_value);
   if minus_log_phi < max_minus_log_phi + increment_minus_log_phi
        minus_log_phi_index = ceil((minus_log_phi - min_minus_log_phi)/increment_minus_log_phi) + 1;
        minus_log_phi_inv_table(minus_log_phi_index) = x_value;
   end
end


x_test = 0:0.01:400;
phi_test2 = phi_x(x_test);
phi_inv_test = phi_x_inv(phi_test2);


figure;
% plot(x_vec, phi_vec, 'LineWidth', 2);
hold on; grid on;
plot(x_test, phi_test2, 'LineWidth', 2);
plot(x_test, phi_inv_test, 'LineWidth', 2);

hold on; grid on;
% hold on; grid on;
% 
% for i = 0:0.01:20
%    plot(i, phi(i), 'r*'); 
% end