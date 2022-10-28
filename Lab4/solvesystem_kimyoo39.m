function [t,x] = solvesystem_kimyoo39(func1,func2,t0,tN,x0,h)
number = round((tN-t0) - h);
tvals = linspace(t0, tN, number);
xvals = zeros(2, number);

for i = 1:2
    xvals(i,1) = x0;
    if i == 1
        for j = 1:number - 1
            left = func1(tvals(j), xvals(i,j));
            right = func1(tvals(j) + h, xvals(i,j) + (h * left));
            xvals(i,j+1) = xvals(i,j) + (h/2) * (left + right);
        end
    elseif i == 2
        for j = 1:number - 1
            left = func2(tvals(i), xvals(i,j));
            right = func2(tvals(i) + h, xvals(i,j) + (h * left));
            xvals(i,j+1) = xvals(i,j) + (h/2) * (left + right);
        end
    end
end
t = tvals;
x = xvals;
end

