%% Header
% This file serves as workbench for the project when workign with the
% lower-dimensional Stockmeyer model. It should also serve as a reasonable
% template for interpreting the project's design.

%% Boilerplate

addpath(genpath(pwd));
parameters = boiler_plate();

disp('Ran Boilerplate')

%% Set up problem

parameters.N_x      = 2^7;
parameters.N_y      = 2^7;
parameters.dt       = 0.001;
parameters.tfin     = 5;
parameters.fr       = 50; 
parameters.pr       = 0.002;
parameters.m_sz     = 15*(2^10)^3/(8*parameters.N_x*parameters.N_y);
parameters.update   = true;
parameters.ef_angle = pi/4;

disp('Set Parameters')

%% Set Initial Conditions

% Set up discretization
x = linspace(-pi, pi, parameters.N_x);
y = linspace(-pi, pi, parameters.N_y);
[X, Y] = meshgrid(x,y);

% Construct initial conditions
y_side = exp(-10*(y+0*pi).^2) + ...
         exp(-10*(y-0.2*pi).^2);
x_side = exp(-10*(x+0*pi).^2);
f = x_side .* y_side';

ic = f/trapz(y, trapz(x, f, 2));

disp('Set Initial Conditions')

%% Visualize initial conditions

initial_frame = figure(1);
clf(initial_frame);
ax = axes(initial_frame);

frame(X,Y,ic,parameters,ax,...
      "Title","Initial State");

disp('Done Visualizing IC')

%% Set interaction

function return_data = A(z,e,f,c)
% Be very careful with the output, it must be of the form [1 ~]

    % No interaction
    %return_data = 0*z;

    % Electric field
    % I changed angle_to_vec while working on 3d problem so this is
    % problaby broken now
    %vec = sum((0-1*e) .* angle_to_perp(z'),2);
    %return_data = vec';

    % Convolution
    % This one works a little differently, rather than outputting the value
    % at a specific angle, it outputs at all angles
    no_vel = trapz(f,1);
    real_interaction = conv(c,no_vel,"same");
    return_data = -real_interaction/max(real_interaction);

end

function return_data = B(z)
    return_data = -1*z;
end

ef = angle_to_vec(parameters.ef_angle);

parameters.A = @(z,ef,f,c) A(z,ef,f,c);
parameters.B = @(z) B(z);

disp('Done defining interaction')

%% Run simulation

[t_c, d_c] = cont_sudospec_etd_trotter(ic,parameters);

disp('Done Simulating')

%% Visualize Final Conditions

final_frame = figure(2);
clf(final_frame);
ax = axes(final_frame);

frame(X,Y,squeeze(d_c(end,:,:)),parameters,ax,...
      "Title","Final State");

%hold on
%m = plot3(squeeze(traj(1,:)), squeeze(traj(2,:)), squeeze(traj(3,:)), 'r');

disp('Done visualizing FC')

%% Animate

animation_frame = figure(3);
clf(animation_frame);
ax = axes(animation_frame);

animate(x,y,d_c,parameters,ax,...
       "Title","Experiment", ...
       "Time",t_c);

disp('Done animation')

%% Test make_video

make_video(x,y,d_c,parameters,figure(4), ...
    'Title',"Stockmayer with Interactive Forcing", ...
    'Time',t_c);

disp('Done making video')