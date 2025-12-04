%% Boilerplate

addpath(genpath(pwd));
parameters = boiler_plate();

disp('Ran Boilerplate')

%% Set up problem

parameters.N_x    = 2^9;
parameters.N_y    = 2^9;
parameters.dt     = 0.0004;
parameters.tfin   = 1;
parameters.fr     = 50; 
parameters.pr     = 0.002;
parameters.m_sz   = 2^13;
parameters.update = true;

disp('Set Parameters')

%% Run tests

%test_horizontal_pipe(parameters);
test_vertical_pipe(parameters);

disp('Ran Tests')

function return_data = test_vertical_pipe(parameters)
    
    x = linspace(-1, 1, parameters.N_x);
    y = linspace(-1, 1, parameters.N_y);
    [X, Y] = meshgrid(x,y);

    x_side = exp(-10*x.^2);
    y_side = 1 + 0*x_side;
    f = x_side .* y_side';
    ic = f/trapz(y, trapz(x, f, 2));

    parameters.A = @(z) -1*z;
    parameters.B = @(z) 0*z;

    [t_c, d_c] = cont_sudospec_etd_trotter(ic,parameters);

    initial_frame = figure(1);
    clf(initial_frame);
    ax = axes(initial_frame);

    frame(X,Y,squeeze(d_c(1,:,:)),parameters,ax,...
          "Title","Initial State");

    no_shear_frame = figure(2);
    clf(no_shear_frame);
    ax = axes(no_shear_frame);

    frame(X,Y,squeeze(d_c(end,:,:)),parameters,ax,...
          "Title","Should Not Shear");

    parameters.A = @(z) 0*z;
    parameters.B = @(z) -1*z;

    [t_c, d_c] = cont_sudospec_etd_trotter(ic,parameters);

    shear_frame = figure(3);
    clf(shear_frame);
    ax = axes(shear_frame);

    frame(X,Y,squeeze(d_c(end,:,:)),parameters,ax,...
          "Title","Should Shear");
    
end

function return_data = test_horizontal_pipe(parameters)
    
    x = linspace(-1, 1, parameters.N_x);
    y = linspace(-1, 1, parameters.N_y);
    [X, Y] = meshgrid(x,y);

    y_side = exp(-10*x.^2);
    x_side = 1 + 0*y_side;
    f = x_side .* y_side';
    ic = f/trapz(y, trapz(x, f, 2));

    parameters.A = @(z) -1*z;
    parameters.B = @(z) 0*z;

    [t_c, d_c] = cont_sudospec_etd_trotter(ic,parameters);

    initial_frame = figure(1);
    clf(initial_frame);
    ax = axes(initial_frame);

    frame(X,Y,squeeze(d_c(1,:,:)),parameters,ax,...
          "Title","Initial State");

    no_shear_frame = figure(2);
    clf(no_shear_frame);
    ax = axes(no_shear_frame);

    frame(X,Y,squeeze(d_c(end,:,:)),parameters,ax,...
          "Title","Should Shear");

    parameters.A = @(z) 0*z;
    parameters.B = @(z) -1*z;

    [t_c, d_c] = cont_sudospec_etd_trotter(ic,parameters);

    shear_frame = figure(3);
    clf(shear_frame);
    ax = axes(shear_frame);

    frame(X,Y,squeeze(d_c(end,:,:)),parameters,ax,...
          "Title","Should Not Shear");
    
end