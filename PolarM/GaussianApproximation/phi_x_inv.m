function [ x ] = phi_x_inv( y )

global minus_log_phi_inv_table;
global min_minus_log_phi;
global max_minus_log_phi;
global increment_minus_log_phi;
minus_log_phi = -log(y);
minus_log_phi = max(minus_log_phi, min_minus_log_phi);
minus_log_phi = min(minus_log_phi, max_minus_log_phi);
minus_log_phi_index = round((minus_log_phi - min_minus_log_phi)/increment_minus_log_phi - 0.499) + 1;
x = minus_log_phi_inv_table(minus_log_phi_index);

end

