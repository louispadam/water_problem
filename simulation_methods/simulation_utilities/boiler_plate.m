function return_data = boiler_plate()

    % Generate parameter structure
    parameters = struct(...
        'N_x',1,...        % x discretization for 2d
        'N_y',1,...        % y discretization for 2d and 3d (bad name)
        'N_omega',1,...    % omega discretization for 3d
        'N_theta',1,...    % theta discretization for 3d
        'A',@(z) 1,...     %
        'B',@(z) 1,...     %
        'del',1,...        % parameter for japanese bracket
        'dt',0.001, ...    % simulation time step
        'tfin',100, ...    % ending time
        'fr',2, ...        % frame rate
        'pr',1, ...        % pause rate
        'm_sz',1, ...      % max vector size
        'update',false, ...   % print progress of simulation
        'ef_angle',0, ...
        'ker_pow',1, ...
        'ef_pow',1);

    return_data = parameters;

end