function return_data = animate_contour(y,omega,theta,S,parameters,fig,options)
%ANIMATE_CONTOUR Animate a countour plot in the high-dimensional model.
%Automatically saves the animation in an mp4 file named output.
%
%last updated 03/10/26 by Adam Petrucci
arguments (Input)
    y                       % discretization in y-coordinate
    omega                   % discretization in omega-coordinate
    theta                   % discretization in theta-coordinate    
    S                       % surface data to plot
    parameters struct       % parameters used for simulation
    fig                     % fig to build in
end
arguments (Input)
    options.Title = ""      % title of axis
    options.Time = []       % time data
    options.xLabel = 'X'    % x-label for figure
    options.yLabel = 'Y'    % y-label for figure
    options.zLabel = 'Z'    % z-label for figure
end

    %****************************
    % Collect Inputs
    %****************************

    clf(fig);
    ax = axes(fig);
    tit = options.Title;
    time = options.Time;

    %****************************
    % Define Temporal Parameters
    %****************************
    dt=parameters.dt;       % simulation time step
    tt=0;                   % current time
    ptfac=parameters.fr;    % frame rate

    %****************************
    % Run Animation
    %****************************

    v = VideoWriter('output.mp4','MPEG-4');
    v.FrameRate = 30;         % Adjust for smoothness
    open(v);

    % Set up annotation for time-keeping
    if ~isempty(time)
        p = ax.Position;
        a = annotation('textbox', ...
            [0.7,0.9, 0.1, 0.1], ...
            'String', sprintf('Time: %d',time(1)), ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', ...
            'FontSize', 11, ...
            'Color','w');
    end

    data_slice = squeeze(S(1,:,:,:));
    cut = 0.9 * max(data_slice,[],"all");

    % There's some meshgrid vs ndgrid ordering nonsense here
    p = patch(ax,isosurface(omega,y,theta,data_slice,cut));
    isonormals(omega,y,theta,data_slice,p);

    set(p,'FaceColor','red','EdgeColor','none');

    daspect(ax,[1 1 1]);
    view(ax,3);
    camlight(ax);
    lighting(ax,'gouraud')

    xlim(ax,[-parameters.w_bound parameters.w_bound])
    ylim(ax,[-parameters.y_bound parameters.y_bound])
    zlim([-pi pi])

    % Parameters for plot
    ylabel(ax,options.yLabel);
    xlabel(ax,options.xLabel);
    zlabel(ax,options.zLabel);
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')

    for ind = 1:length(time)

        if mod(ind, ptfac) == 0 % I should be able to speed this up by
                                % putting it in the for loop

            % Collect slices of data for current frame
            data_slice = squeeze(S(ind,:,:,:));
            cut = 0.9 * max(data_slice,[],"all");

            % Display frame for given time
            p = patch(ax,isosurface(omega,y,theta,data_slice,cut));
            isonormals(omega,y,theta,data_slice,p)
            set(p,'FaceColor','red','EdgeColor','none')

            % Update annotation tracking time
            if ~isempty(time)
                a.String = sprintf('Time: %.2f',time(ind));
            end
            
            writeVideo(v, getframe(fig));

        end
        tt=tt+dt; % Update time
    end

    return_data = 1;
end