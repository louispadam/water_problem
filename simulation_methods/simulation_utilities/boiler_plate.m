function return_data = boiler_plate()

    % Generate parameter structure
    parameters = struct(...
        'N_x',1,...        % x discretization
        'N_y',1,...        % y discretization
        'A',@(z) 1,...     %
        'B',@(z) 1,...     %
        'dt',0.001, ...    % simulation time step
        'tfin',100, ...    % ending time
        'fr',2, ...        % frame rate
        'pr',1, ...        % pause rate
        'm_sz',1, ...      % max vector size
        'update',false, ...   % print progress of simulation
        'ef_angle',0);

    return_data = parameters;

end