function return_data = frame(X,Y,S,parameters,axis,options)
%FRAME Present simulation data on a provided axis handle. Accepts a cell
% array, each element of which is a vector of particles.
%
%last updated 10/07/25
arguments (Input)
    X                       % mesh for x-coordinate
    Y                       % mesh for y-coordinate
    S                       % surface to plot
    parameters struct       % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Title = ""                      % title of axis
end

    %****************************
    % Collect Inputs
    %****************************

    ax = axis;
    tit = options.Title;

    %****************************
    % Construct Figures
    %****************************

    % Plot
    surf(X, Y, S);
    shading interp;
    view(2);
    colormap(parula);
    colorbar;

    % Parameters for plot
    ylabel(ax,'Y');
    xlabel(ax,'X');
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')

    return_data = ax;

end