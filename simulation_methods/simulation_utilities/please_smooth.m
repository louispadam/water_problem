function return_data = please_smooth(x,y,f,delta)
%PLEASE_SMOOTH Smooth out an input function f on an (x)x(y) grid using
%convolution against a Gaussian kernel.
%
%last updated 03/09/26 by Adam Petrucci

    % Enlarge f to account for periodic boundaries
    % This could be made more efficient: extend by some function of the
    % smoothing parameter delta rather than the full entire space
    f_p = repmat(f,3,3);

    % Constructing smoothing kernel
    ex = @(z1,z2) exp(-(z1^2+z2^2)/delta);
    e = arrayfun(ex,x,y);

    % Convolve function with kernel
    smoothed = conv2(f_p,e,'same');

    % Extract original dimensions
    s = size(x);
    smoothed_og = smoothed(s(1):2*s(1)-1,s(2):2*s(2)-1);

    return_data = smoothed_og;

end