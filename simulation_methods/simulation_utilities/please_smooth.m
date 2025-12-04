function return_data = please_smooth(X,Y,f,delta)

    % Enlarge f to account for periodic boundaries
    % This could be made more efficient if I incorporated the smoothing
    % parameter delta
    f_p = repmat(f,3,3);

    % Constructing smoothing kernel
    ex = @(z1,z2) exp(-(z1^2+z2^2)/delta);
    e = arrayfun(ex,X,Y);

    % Convolve function with kernel
    smoothed = conv2(f_p,e,'same');

    % Extract original dimensions
    s = size(X);
    smoothed_og = smoothed(s(1):2*s(1)-1,s(2):2*s(2)-1);

    return_data = smoothed_og;

end