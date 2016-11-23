
function [y] = phi_x(x)

y = (x <  10) .* exp(-0.4527 * x.^(0.86) + 0.0218);
y = (x >= 10) .* sqrt(pi./x) .* (1 - 1.4286./(x + 0.0001)) .* exp(-x/4) + y;
% y = min(1, y);