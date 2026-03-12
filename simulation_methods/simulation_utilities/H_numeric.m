function return_data = H_numeric(params)
%H_NUMERIC Numerical calculation of the convolution kernel. In particular,
%the exact solution integrates out a spatial dimension exactly, whereas
%it's computed numerically here.
%
%last updated 03/11/26 by Adam Petrucci

    L = 2*params.y_bound;

    % set up important discretizations
    % discretization in x-dimension is significantly finer than in y to
    % reduce error from quadrature
    y = linspace(-params.y_bound,params.y_bound,params.N_y);
    x = linspace(-params.y_bound,params.y_bound,128*params.N_y);

    x_dists = periodic_influence(0,x,L);
    y_dists = periodic_influence(0,y,L);

    [conv_X, conv_Y] = ndgrid(x_dists,y_dists);

    % compute four terms in matrix
    rr2 = (conv_X.^2 + conv_Y.^2 + params.del^2);
    
    front_coeff = 1./rr2.^2;
    top_left = rr2 - 2*conv_X.^2;
    off_diag = -2*conv_X.*conv_Y;
    bottom_right = rr2 - 2*conv_Y.^2;

    ker_tl = front_coeff.*top_left;
    ker_of = front_coeff.*off_diag;
    ker_br = front_coeff.*bottom_right;
   
    % integrate out spatial x-dimension
    out_x_tl = trapz(x,ker_tl,1);
    out_x_of = trapz(x,ker_of,1);
    out_x_br = trapz(x,ker_br,1);

    % construct matrix
    ker = [out_x_tl; out_x_of; out_x_of; out_x_br]';

    return_data = ker;

end