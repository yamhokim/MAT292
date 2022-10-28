function [x,y] = odesolver(func,t0,tN,y0,h)
% specify func for later steps
number = (tN-t0)/h; % Create the number of steps
tvals = linspace(t0,tN,number); % Create an array of equally sized steps from t0 to tN
yvals = zeros(size(tvals)); % Create an array of zeros to be filled with yvals
yvals(1) = y0; % set the first yval in yvals to y0
for i = 1:(number-1) % Run a for loop to compute the yval at each step using improved eulers method
    left = func(tvals(i), (yvals(i)));
    right = func(tvals(i)+h, yvals(i)+(h*left));
    yvals(i+1) = yvals(i) + (h/2)*(left + right);
end
x = tvals; % output the tvals/xvals
y = yvals; % output corresponding yvals
end
