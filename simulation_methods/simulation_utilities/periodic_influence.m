function return_data = periodic_influence(x,y,L)
%PERIODIC_INFLUENCE Calculates periodic distance between inputs x and y on
%an interval of length L. Designed to accept vectors and return the
%corresponding grid of distances.
%
%last updated 03/09/26 by Adam Petrucci

    orientation = -sign(x-y') .* sign(L/2 - abs(x-y'));
    magnitude = (L/2 - abs(L/2 - abs(x-y')));

    % If points are maximally separated then orientation is ambiguous and
    % must be manually set to recognize any distance.
    orientation(orientation == 0) = 1;

    return_data = orientation .* magnitude;

end