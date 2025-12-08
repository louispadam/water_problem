function return_data = periodic_influence(x,y,L)
    orientation = -sign(x-y') .* sign(L/2 - abs(x-y'));
    magnitude = (L/2 - abs(L/2 - abs(x-y')));
    return_data = orientation .* magnitude;
end