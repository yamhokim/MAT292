function [x,y] = odesolver(func,t0,tN,y0,h)
% specify func for later steps
number = tN/h;
tvals = linspace(t0,tN,number);
yvals = zeros(size(tvals));
yvals(1) = y0;
for i = 1:(number-1)
    left = func(tvals(i), (yvals(i)));
    right = func(tvals(i+1), yvals(i));
    yvals(i+1) = yvals(i) + (h/2)*(left + right);
end
x = tvals;
y = yvals;
end