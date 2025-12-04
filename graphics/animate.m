function return_data = animate(X,Y,S,parameters,axis,options)
%ANIMATE Animate a collection of simulations. Accepts a cell array, each
%element of which is an 2d array of data: time x particle.
%
%last updated 10/07/25 by Adam Petrucci
arguments (Input)
    X                       % mesh for x-coordinate
    Y                       % mesh for y-coordinate
    S                       % surface to plot
    parameters struct       % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Title = ""      % title of axis
    options.Time = []
end

    %****************************
    % Collect Inputs
    %****************************

    ax = axis;
    tit = options.Title;
    time = options.Time;

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

    for ind = 1:length(time)

        if mod(ind, ptfac) == 0 % I should be able to speed this up by
                                % putting it in the for loop

            cla(ax,'reset');

            % Collect slices of data for current frame
            to_send = squeeze(S(ind,:,:));

            % Display frame for given time
            frame(X,Y,to_send,parameters,ax, ...
                  Title = tit);

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