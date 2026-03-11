function return_data = contour(y,theta,omega,S,parameters,axis,options)
%CONTOUR Constructs a contour plot in the high-dimensional model.
%
%last updated 03/10/26 by Adam Petrucci
arguments (Input)
    y                       % mesh for y-coordinate
    theta                   % mesh for theta-coordinate
    omega                   % mesh for omega-coordinate
    S                       % surface data to plot
    parameters struct       % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Title = ""      % title of axis
    options.xaxis = "X"     % x-label of figure
    options.yaxis = "Y"     % y-label of figure
    options.zaxis = "Z"     % z-label of figure
end

    %****************************
    % Collect Inputs
    %****************************

    ax = axis;
    tit = options.Title;
    xlab = options.xaxis;
    ylab = options.yaxis;
    zlab = options.zaxis;

    %****************************
    % Construct Figures
    %****************************

    cut = 0.9 * max(S,[],"all");

    % There's some meshgrid vs ndgrid ordering nonsense here
    p = patch(ax,isosurface(omega,y,theta,S,cut));
    isonormals(omega,y,theta,S,p)

    set(p,'FaceColor','red','EdgeColor','none')

    daspect(ax,[1 1 1])
    view(ax,3)
    camlight(ax)
    lighting(ax,'gouraud')

    xlim(ax,[-pi pi])
    ylim(ax,[-pi pi])
    zlim(ax,[-pi pi])

    % Parameters for plot
    if ylab~="none"
        ylabel(ax,ylab);
    end
    if xlab~="none"
        xlabel(ax,xlab);
    end
    if zlab~="none"
        zlabel(ax,zlab);
    end
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')
    xlim(ax,[-pi pi]);
    ylim(ax,[-pi pi]);

    return_data = ax;

end