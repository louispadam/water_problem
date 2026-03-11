function return_data=track_peak(x,y,data,parameters)
%TRACK_PEAK In the case of two spatial dimensions, track the point of
%maximum intensity.
%
%last updated 03/09/26 by Adam Petrucci
arguments
    x           % discretization in x dim
    y           % discretization in y dim
    data        % time series of simulation.
    parameters  % parameters for experiment
end

    s = size(data);

    % Find maximum signal across all time steps
    [val, idx] = max(reshape(data,s(1),[]), [], 2);
    [row, col] = ind2sub([parameters.N_x parameters.N_y], idx);

    % Construct artificial z-value so that trajectory will plot above data
    % if presented in the same figure.
    trak_max_x = x(col);
    trak_max_y = y(row);
    z_placeholder = max(val)*ones([1,s(1)])+1;
    trajectory = cat(1,trak_max_x,cat(1,trak_max_y,z_placeholder));

    return_data = trajectory;

end
