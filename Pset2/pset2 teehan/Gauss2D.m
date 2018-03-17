% x0 = [x y sigma amp b]


function strf = Gauss2D(x0, xdata)
strf = x0(5) + x0(4) * ...
    exp(-(((xdata(:,1) - x0(1)).^2 + (xdata(:,2) - x0(2)).^2) / (2 * x0(3)^2)));
end