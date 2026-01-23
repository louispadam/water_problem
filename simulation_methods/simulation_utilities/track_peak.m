function return_data=track_peak(x,y,data,parameters)
%
%last updated 01/19/26 by Adam Petrucci
arguments
    x
    y
    data
    parameters
end

    s = size(data);

    [val, idx] = max(reshape(data,s(1),[]), [], 2);
    [row, col] = ind2sub([parameters.N_x parameters.N_y], idx);

    trak_max_x = x(col);
    trak_max_y = y(row);
    z_placeholder = max(val)*ones([1,s(1)])+1;

    trajectory = cat(1,trak_max_x,cat(1,trak_max_y,z_placeholder));

    return_data = trajectory;

end
