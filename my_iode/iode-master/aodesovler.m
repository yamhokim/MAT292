function [tvals,yvals] = aodesolver(f, t0, tN, y0, h)
    tol = 1e-8;
    yvals = [y0];
    tvals = [t0];
    n = 2;

    while (tvals(length(tvals)) < tN)
        y = f(tvals(n-1), yvals(n-1))*h +yvals(n-1);
        z_1 = f(tvals(n-1), yvals(n-1))*h + yvals(n-1);
        z_2 = f(tvals(n-1)+h/2, z_1)*(h/2) + z_1;
        diff = abs(z_2 - y);
        if (diff < tol)
            yvals = [yvals z_2+(z_2-y)];
            tvals = [tvals tvals(n-1)+h];
            n = n+1;
            break;
        else
            h = 0.9*h*min(max(tol/abs(z_2-y),0.3),2);
        end
    end
end