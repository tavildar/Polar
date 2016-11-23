function [ y ] = phi_x_table( x )

global phi_x_table;
global min_x;
global max_x;
global increment_x;
x = max(min_x, x);
x = min(max_x, x);
x_index = round((x - min_x)/increment_x) + 1;
y = phi_x_table(x_index);

end

