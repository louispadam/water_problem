%% Header
% This file serves as workbench for the project when workign with the full
% Stockmeyer model. It should also serve as a reasonable template for
% interpreting the project's design.

%% Boilerplate

addpath(genpath(pwd));
parameters = boiler_plate();

disp('Ran Boilerplate')

%% Set up problem

parameters.N_y      = 2^5;
parameters.N_omega  = 2^6;
parameters.N_theta  = 2^7;
parameters.dt       = 0.005;
parameters.tfin     = 10;
parameters.fr       = 10; 
parameters.pr       = 0.002;
parameters.del      = 0.001;
parameters.m_sz     = 15*(2^10)^3/(8*parameters.N_y*...
                                     parameters.N_omega*...
                                     parameters.N_theta);
parameters.update   = true;
parameters.ef_angle = pi/2;
parameters.ker_pow  = 1/100;
parameters.ef_pow   = 0;

disp('Set Parameters')

%% Set Initial Conditions

% Set up discretization
y = linspace(-pi, pi, parameters.N_y);
omega = linspace(-pi, pi, parameters.N_omega);
theta = linspace(-pi,pi, parameters.N_theta);

[Omega_2, Theta_2] = ndgrid(omega,theta);
[Y_3, Omega_3, Theta_3] = ndgrid(y,omega,theta);

% Construct initial conditions in each dimension
y_side     = exp(-5*(y+0*pi).^2);
omega_side = exp(-10*(omega+0*pi).^2);
theta_side = exp(-5*(theta+0*pi).^2) + ...
             exp(-5*(theta-0.3*pi).^2);

% Tensor-product them together
f = y_side' .* omega_side .* reshape(theta_side,1,1,[]);

% Normalize initial conditions
ic = f/trapz(y, trapz(omega, trapz(theta, f, 3), 2));

disp('Set Initial Conditions')

%% Visualize initial middle slice

initial_frame = figure(1);
clf(initial_frame);
ax = axes(initial_frame);

frame(Omega_2,Theta_2,squeeze(ic(round(parameters.N_y/2),:,:)),parameters,ax,...
      "Title","Initial State",...
      "xaxis",'\theta',...
      "yaxis",'\omega');

disp('Done Visualizing IC')

%% Set interaction

conv_data = H_exact(parameters);

parameters.A = @(f) A(theta,f,conv_data,parameters);
parameters.B = @(z) B(z);

disp('Done defining interaction')

%% Visualize convolution kernel

parameters.del = 0.001;

parameters.N_y = 2^5;
y_coarse = linspace(-pi, pi, parameters.N_y);
H_data_coarse = H_exact(parameters);

parameters.N_y = 2^9;
y_fine = linspace(-pi, pi, parameters.N_y);
H_data_fine = H_exact(parameters);

convolution_coarse_frame = figure(1);
clf(convolution_coarse_frame);
ax = axes(convolution_coarse_frame);

hold on
plot(ax,y_coarse,H_data_coarse(:,1),...
     "LineWidth",2.0,...
     "DisplayName",'(1,1)-entry (h_{11})')
plot(ax,y_coarse,H_data_coarse(:,4),...
     "LineWidth",2.0,...
     "DisplayName",'(2,2)-entry (h_{22})')

title(ax,'Kernel with Coarse Mesh (N=2^5, \delta=0.001)')
xlabel(ax,'y (spatial dimension)')
legend(ax)

convolution_fine_frame = figure(2);
clf(convolution_fine_frame);
ax = axes(convolution_fine_frame);

hold on
plot(ax,y_fine,H_data_fine(:,1),...
     "LineWidth",2.0,...
     "DisplayName",'(1,1) entry (h_{11})')
plot(ax,y_fine,H_data_fine(:,4),...
     "LineWidth",2.0,...
     "DisplayName",'(2,2) entry (h_{22})')

title(ax,'Kernel with Fine Mesh (N=2^9, \delta=0.001)')
xlabel(ax,'y (spatial dimension)')
legend(ax)

saveas(figure(1),'convolution_coarse_mesh.png')
saveas(figure(2),'convolution_fine_mesh.png')

disp('Plotted convolution kernels')

%% Visualize continous and singular components

function return_data = get_continuous_piece(params)

    L = 2*pi;
    d = params.del;
    fi = linspace(-pi,pi,params.N_y);
    inp = periodic_influence(0,fi,2*pi);
    jb = jap_brac(inp,d);
    t1f1 = (d^2-inp.^2)./(jb.^2) - 1;
    t1f2 = (2*L./(L^2+4*jb.^2));

    return_data = t1f1.*t1f2;

end

function return_data = get_singular_piece(params)

    L = 2*pi;
    d = params.del;
    fi = linspace(-pi,pi,params.N_y);
    inp = periodic_influence(0,fi,2*pi);
    jb = jap_brac(inp,d);
    t2f1 = 1+((d^2-inp.^2)./(jb.^2));
    t2f2 = atan(L./(2*jb))./jb;

    return_data = t2f1.*t2f2;

end

parameters.N_y = 2^15;

cont_piece = get_continuous_piece(parameters);
sing_piece = get_singular_piece(parameters);


convolution_cont_frame = figure(1);
clf(convolution_cont_frame);
ax = axes(convolution_cont_frame);

hold on
plot(linspace(-pi,pi,parameters.N_y),cont_piece,...
     "LineWidth",2.0)

title(ax,'First component of kernel')
xlabel(ax,'y (spatial dimension)')

convolution_sing_frame = figure(2);
clf(convolution_sing_frame);
ax = axes(convolution_sing_frame);

hold on
plot(linspace(-pi,pi,parameters.N_y),sing_piece,...
     "LineWidth",2.0)

title(ax,'Second component of kernel')
xlabel(ax,'y (spatial dimension)')

saveas(figure(1),'convolution_first_component.png')
saveas(figure(2),'convolution_second_component.png')

disp('Done splitting convolution kernel')

%% Run simulation

[t_c, d_c] = cont_sudospec_etd_trotter_3d(ic,parameters);

disp('Done Simulating')

%% Animate contour plot

animate_contour(y,omega,theta,d_c,parameters,figure(4), ...
    'Title',sprintf("testing testing"), ...
    'Time',t_c,...
    'xLabel','\theta',...
    'yLabel','\omega',...
    'zLabel','another');

disp('Done making contour video')

%% Visualize final middle slice

initial_frame = figure(5);
clf(initial_frame);
ax = axes(initial_frame);

frame(Omega_2,Theta_2,squeeze(d_c(end,round(parameters.N_y/2),:,:)),parameters,ax,...
      "Title","Final State (everything)", ...
      'xaxis','\theta',...
      'yaxis','\omega');

saveas(figure(5),'final.png')

disp('Done Visualizing FC')

%% Make video of middle slice

make_video(omega,theta,squeeze(d_c(:,round(parameters.N_y/2),:,:)),parameters,figure(4), ...
    'Title','Test Without Electric Field (\nu=0.01, \delta=0.001)', ...
    'Time',t_c,...
    'xLabel','\theta',...
    'yLabel','\omega');

disp('Done making middle-slice video')
