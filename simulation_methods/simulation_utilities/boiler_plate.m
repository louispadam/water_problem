function return_data = boiler_plate()
%BOILER_PLATE Generates a structure containing all parameters for a
%simulation of water. I recommend calling this function at the beginning of
%any experiment and modifying parameter values from there.
%
%last updated 03/11/26 by Adam Petrucci

    % Generate parameter structure
    parameters = struct(...
        'y_bound',...         % bound in spatial y-dim
        'w_bound',...         % bound in omega
        'N_x',1,...           % x discretization in 1 spatial dim
        'N_y',1,...           % y discretization in 2 and 3 spatial dim
        'N_omega',1,...       % omega discretization in 3 spatial dim
        'N_theta',1,...       % theta discretization in 3 spatial dim
        'A',@(z) 1,...        % theta coefficient of pde
        'B',@(z) 1,...        % omega coefficient of pde
        'del',1,...           % parameter for japanese bracket
        'dt',0.001, ...       % simulation time step
        'tfin',100, ...       % ending time
        'fr',2, ...           % frame rate (for animation)
        'pr',1, ...           % pause rate (for animation)
        'm_sz',1, ...         % max vector size (for data storage)
        'update',false, ...   % print progress of simulation
        'ef_angle',0, ...     % angle of electric field
        'ker_pow',1, ...      % coefficient on interaction term of pde
        'ef_pow',1);          % coefficient on electric field term of pde

    return_data = parameters;

end