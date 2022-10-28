function [tvals,yvals] = aodesolver(f, t0, tN, y0, h)
    bound = 1e-8; % Bound for error
    tvals = [t0]; % Create an array initially only containing t0, this will hold all future tvals
    yvals = [y0]; % Create an array initially only containing y0, this will hold all future yvals
    n = 2; % used for indexing

    while (tvals(length(tvals)) < tN) % While the newest tval is less than the upper bound of tvals, run the while loop
        while 1 % Infinitely running while loop that is only stopped when our diff < bound
            y = f(tvals(n-1), yvals(n-1))*h + yvals(n-1); % Calculate the yval
            z_1 = f(tvals(n-1), yvals(n-1))*(h/2) + yvals(n-1);
            z_2 = f(tvals(n-1)+(h/2), z_1)*(h/2) + z_1; % calculate the second "yval"
            diff = abs(z_2-y); % determine the difference between the two values
            if (diff < bound) % The error is within the tolerance range, so no change needs to be made to stepsize. 
                yvals(end+1) = z_2 + (z_2-y); % concatenate the latest yval to the end of the yvals array
                tvals(end+1) = tvals(n-1)+h; % concatenate the latest tval to the end of the tvals array
                n = n+1; % increase the index and move onto the next step
                break
            elseif (diff >= bound)% The error is larger than the tolerance range, so the stepsize must be changed.
                h = 0.9*h*min(max(bound/abs(z_2-y),0.3),2); % update the stepsize and recalculate the y, z_1, z_2, and diff, then reevaluate if the stepsize is valid
            end
        end
    end
