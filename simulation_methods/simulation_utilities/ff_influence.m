function return_data = ff_influence(x)
    shape = exp(-5*(x+pi/16).^2) - exp(-5*(x-pi/16).^2);
    %shape = sech(3*x);
    return_data = shape/max(shape);
end