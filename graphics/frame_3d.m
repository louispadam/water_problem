function return_data = frame_3d(fig_num, data, params, options)
%FRAME_3d Present simulation data in omega, theta, and y by fitting regular
%(low dimensional) frames into a tiled layout, each tile representing a y
%value.
%
%last updated 03/10/26
arguments (Input)
    fig_num                    % figure number
    data                       % data to plot
    params                     % experiment parameters
end
arguments (Input)
    options.Title  = "3D Visualization"    % collection title
    options.xaxis  = "Theta"               % x-axis label
    options.yaxis  = "Omega"               % y-axis label
    options.layers = 12;                   % number of slices
end

    layers = options.layers;

    fig = figure(fig_num);
    fig.Units = 'normalized';
    fig.Position = [0.01 0.1 0.98 0.75];

    tl = tiledlayout("flow");
    tl.Padding = 'compact';
    tl.TileSpacing = 'compact';

    omega = linspace(-params.w_bound, params.w_bound, params.N_omega);
    theta = linspace(-pi,pi, params.N_theta);
    [Omega_2, Theta_2] = ndgrid(omega,theta);

    M = max(data,[],'all');

    % cycle through frames
    for i = 1:layers
        ax = nexttile;
        layer = floor((i-1)*params.N_y/layers)+1;
        frame(Omega_2,Theta_2,squeeze(data(layer,:,:)),params,ax,...
        "Title",sprintf('Slice %0.2f',i/layers),...
        "xaxis","none",...
        "yaxis","none");
        ax.CLim = [0,M];
    end

    title(tl,options.Title)
    xlabel(tl,options.xaxis)
    ylabel(tl,options.yaxis)

    return_data = 1;

end