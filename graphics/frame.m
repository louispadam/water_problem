function return_data = frame(X,Y,S,parameters,axis,options)
%FRAME Present simulation data in only omega and theta on a provided axis
%handle.
%
%last updated 03/10/26
arguments (Input)
    X                       % mesh for x-coordinate
    Y                       % mesh for y-coordinate
    S                       % surface to plot
    parameters struct       % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Title = ""      % title of axis
    options.xaxis = "X"     % x-label of figure
    options.yaxis = "Y"     % y-label of figure
end

    %****************************
    % Collect Inputs
    %****************************

    ax = axis;
    tit = options.Title;
    xlab = options.xaxis;
    ylab = options.yaxis;

    %****************************
    % Construct Figures
    %****************************

    % Plot
    surf(Y, X, S);
    shading interp;
    view(2);
    colormap(parula);
    colorbar;

    % Parameters for plot
    if ylab~="none"
        ylabel(ax,ylab);
    end
    if xlab~="none"
        xlabel(ax,xlab);
    end
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')
    % maybe be problems here if we return to lower dim simulation
    xlim(ax,[-pi pi]);
    ylim(ax,[-parameters.w_bound parameters.w_bound]);

    return_data = ax;

end