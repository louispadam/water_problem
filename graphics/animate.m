function return_data = animate(x,y,S,parameters,axis,options)
%ANIMATE Animate data in omega and theta (so no spatial representation).
%Everything gets plotted into inputted axis.
%
%last updated 03/10/26 by Adam Petrucci
arguments (Input)
    x                       % discretization in x-coordinate
    y                       % discretization in y-coordinate
    S                       % surface to plot
    parameters struct       % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Title = ""      % title of axis
    options.Time = []       % time data
    options.Trajectory = [] % also plot a given trajectory
end

    %****************************
    % Collect Inputs
    %****************************

    ax = axis;
    tit = options.Title;
    time = options.Time;
    tr = options.Trajectory;

    %****************************
    % Define Temporal Parameters
    %****************************
    dt=parameters.dt;       % simulation time step
    tt=0;               % current time
    ptfac=parameters.fr;    % frame rate
    pr = parameters.pr;     % pause rate

    %****************************
    % Run Animation
    %****************************

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
    ylim(ax,[-parameters.w_bound parameters.w_bound]);
    ylabel(ax,'Y');
    xlabel(ax,'X');
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')

    for ind = 1:length(time)

        if mod(ind, ptfac) == 0 % I should be able to speed this up by
                                % putting it in the for loop

            % Collect slices of data for current frame
            to_send = squeeze(S(ind,:,:));

            % Display frame for given time
            s.ZData = to_send;

            if ~isempty(tr)
                m.XData = tr(1,1:ind);
                m.YData = tr(2,1:ind);
                m.ZData = tr(3,1:ind);
            end

            drawnow limitrate

            % Update annotation tracking time
            if ~isempty(time)
                a.String = sprintf('Time: %.2f',time(ind));
            end
            
            pause(pr);

        end
        tt=tt+dt; % Update time
    end

    return_data = 1;
end