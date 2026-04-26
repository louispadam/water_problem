function return_data = H_lim(params)
%H1_LIM Convolution data in the delta->0 limit excluding point mass
%
%last updated 03/16/26 by Adam Petrucci

    y = linspace(-params.y_bound,params.y_bound,params.N_y);
    r = periodic_influence(0,y,2*pi);
    
    L = 2*params.y_bound;

    num = 4*L;
    denom = L^2+4*r.^2;
    entry = (num./denom);

    return_data = [entry,zeros(size(r)),...
                   zeros(size(r)),-entry];

end