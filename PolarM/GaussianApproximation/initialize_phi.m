function  initialize_phi(x_increment)

persistent phi_initialized;
global minus_log_phi_inv_table;
global min_minus_log_phi;
global max_minus_log_phi;
global increment_minus_log_phi;
global phi_x_table;
global min_x;
global max_x;
global increment_x;

if isempty(phi_initialized)
    
    phi_initialized = 1;
    if nargin < 1
        x_increment = 10^(-4);
    end
    x_vec = 0: x_increment: 400;
    
    
    phi_vec = (x_vec <  10) .* exp(-0.4527 * x_vec.^(0.86) + 0.0218);
    phi_vec = (x_vec >= 10) .* sqrt(pi./(x_vec+0.0001)) .* (1 - 1.4286./(x_vec + 0.0001)) .* exp(-x_vec/4) + phi_vec;
    phi_vec = min(phi_vec, 1);
    min_minus_log_phi = 0;
    max_minus_log_phi = 100;
    increment_minus_log_phi = 10^(-3);
    minus_log_phi_inv_table = zeros(ceil((max_minus_log_phi - min_minus_log_phi)/increment_minus_log_phi) + 1, 1);
    
    min_x = 0;
    max_x = 100;
    increment_x = 0.01;
    index = 1;
    for x = min_x: increment_x: max_x + increment_x
        if (x < 10)
            phi_x_table(index) = exp(-0.4527 * x.^(0.86) + 0.0218);
        else
            phi_x_table(index) =sqrt(pi./(x)) .* (1 - 1.4286./(x)) .* exp(-x/4);
        end
        index = index + 1;
    end
    
    
    for x_index = 1 : length(x_vec)
        x_value  = x_vec(x_index);
        phi_value = phi_vec(x_index);
        minus_log_phi = -log(phi_value);
        if minus_log_phi < max_minus_log_phi + increment_minus_log_phi
            minus_log_phi_index = ceil((minus_log_phi - min_minus_log_phi)/increment_minus_log_phi) + 1;
            minus_log_phi_inv_table(minus_log_phi_index) = x_value;
        end
    end
    
end

end

