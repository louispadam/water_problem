function return_data = ff_influence(x)
%FF_INFLUENCE Given a discrectized interval, compute an interaction kernel.
%This was used to generate the kernel for experiments in one spatial
%dimension.
%
%last updated 03/09/26 by Adam Petrucci

    shape = exp(-5*(x+pi/16).^2) - exp(-5*(x-pi/16).^2);
    return_data = shape/max(shape);

end