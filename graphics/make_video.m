function return_data = make_video(x,y,S,parameters,fig,options)
%ANIMATE Animate a collection of simulations. Accepts a cell array, each
%element of which is an 2d array of data: time x particle.
%
%last updated 10/07/25 by Adam Petrucci
arguments (Input)
    x                       % discretization in x-coordinate
    y                       % discretization in y-coordinate
    S                       % surface to plot
    parameters struct       % parameters used for simulation
    fig                     % figure to work from
end
arguments (Input)
    options.Title = ""      % title of axis
    options.Time = []
    options.Trajectory = []
end

    %****************************
    % Collect Inputs
    %****************************

    clf(fig);
    ax = axes(fig);
    tit = options.Title;
    time = options.Time;
    tr = options.Trajectory;

    %****************************
    % Define Temporal Parameters
    %****************************
    dt=parameters.dt;       % simulation time step
    tt=0;               % current time
    ptfac=parameters.fr;    % frame rate

    %****************************
    % Run Animation
    %****************************

    v = VideoWriter('output.mp4','MPEG-4');
    v.FrameRate = 30;         % Adjust for smoothness
    open(v);

    [X, Y] = meshgrid(x,y);

    % Set up annotation for time-keeping
    if ~isempty(time)
        p = ax.Position;
        a = annotation('textbox', ...
            [p(1)+0.7*p(3),p(2)+0.9*p(4), 0.1, 0.1], ...
            'String', sprintf('Time: %d',time(1)), ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', ...
            'FontSize', 11, ...
            'Color','w');
    end

    s = surf(X, Y, squeeze(S(1,:,:)));
    shading interp;
    view(2);
    colormap(parula);
    colorbar;

    if ~isempty(tr)
        
        hold on
        m = plot3(tr(1,1), tr(2,1), tr(3,1), 'r');

    end

    % Parameters for plot
    xlim(ax,[-pi pi]);
    ylim(ax,[-pi pi]);
    ylabel(ax,'Y');
    xlabel(ax,'X');
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')

    for ind = 1:length(time)

        if mod(ind, ptfac) == 0 % I should be able to speed this up by
                                % putting it in the for loop

            %cla(ax,'reset');

            % Collect slices of data for current frame
            to_send = squeeze(S(ind,:,:));

            % Display frame for given time
            s.ZData = to_send;

            if ~isempty(tr)
                m.XData = tr(1,1:ind);
                m.YData = tr(2,1:ind);
                m.ZData = tr(3,1:ind);
            end

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