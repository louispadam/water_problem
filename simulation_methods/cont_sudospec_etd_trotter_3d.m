function [return_time, return_data]=cont_sudospec_etd_trotter_3d(initial,parameters)
%
%last updated 11/07/25 by Adam Petrucci
arguments
    initial             % initial conditions
    parameters struct   % parameters for simulation
end

    % Begin timer
    tic

    %****************************
    % Collect Inputs
    %****************************
    ic = initial;

    %****************************
    % System Parameters
    %****************************
    msz = parameters.m_sz;
    ud = parameters.update;
    n_y = parameters.N_y;
    n_omega = parameters.N_omega;
    n_theta = parameters.N_theta;

    %******************************
    % Set up Fourier Transform
    %******************************

    omega = linspace(-pi, pi, n_omega);
    theta = linspace(-pi, pi, n_theta);

    domega = omega(2) - omega(1);
    dtheta = theta(2) - theta(1);

    freq_omega = fftshift( (-n_omega/2 : n_omega/2-1) / (n_omega*domega) );
    freq_theta = fftshift( (-n_theta/2 : n_theta/2-1) / (n_theta*dtheta) );

    %***************************************************
    % Set up iteration
    %**************************************************
    dt=parameters.dt;
    t_final = parameters.tfin;
    U=ic;
    tt=0;

    % Define stepping for iteration
    steps = round(t_final/dt + 1);
    sz = steps;
    keep = 1;   % frequency with which to store iteration

    % If default time-vector is longer than permitted, replace with max
    if sz > msz
        sz = msz;
        keep = steps/msz;
    end

    % Define time and space vectors to store
    si = size(ic);
    time = zeros([1,sz]);
    data = zeros([sz,si]);
    data(1,:,:,:) = ic;
    
    k = 2;      % counter for storing iteration
    here = round(keep*k);
    
    % Prepare frequency of updates (if desired)
    pb = round(linspace(2,steps,20));
    n_pb = 1;

    %***************************************************
    % Iterate!
    %**************************************************

    % updates if desired
    if ud
        fprintf("Began Simulation\n");
    end

    Bomega = parameters.B(omega);
    %temp1 = exp(1i*dt*freq_theta.*Bomega');
    %temp2 = reshape(temp1,1, 64, 128);
    %sg_omega = repmat(temp2,32,1,1);
    %size(sg_omega)
    sg_omega = reshape(exp(1i*dt*freq_theta.*Bomega'),...
                       1, n_omega, n_theta) .* ones(n_y, 1, 1);
    %sg_omega = exp(1i*dt*freq_theta.*Bomega');

    temp = [];

    for step = 2:steps

        % take Fourier in theta
        X_n0 = fft(U,[],3);

        % apply first semigroup (interesting factor in omega)
        X_n12 = X_n0 .* sg_omega;

        % invert Fourier in theta and compute second semigroup
        U = real(ifft(X_n12,[],3));

        Atheta = parameters.A(U);
        temp(end+1) = sum(abs(Atheta),"all");
        sg_theta = exp(1i*dt*reshape(Atheta,n_y,1,n_theta).*...
                             reshape(freq_omega,1,n_omega,1));
        %sg_theta = exp(1i*dt*Atheta.*freq_omega');
        %size(sg_theta)


        % apply Fourier in omega
        Y_n12 = fft(U,[],2);

        % apply second semigroup (interesting factor in theta)
        Y_n1 = Y_n12 .* sg_theta;

        % invert Fourier in omega
        U = ifft(Y_n1,[],2);

        tt=tt+dt;

        % Store result at previously calculated frequency
        if step == here
            time(k) = tt;
            data(k,:,:,:) = real(U);
            k = k+1;
            here = round(keep*k);
        end

        % display update if desired
        if ud && step == pb(n_pb)
            fprintf('%d  Simulation Progress: %3.0f%%\n',step,100*step/steps)
            n_pb = n_pb + 1;
        end

    end

    % Return data
    return_time = time;
    return_data = data;

    % update if desired
    if ud
        fprintf('Completed Simulation in %f seconds\n',toc)
    end

    use_frame = figure(8);
    clf(use_frame);
    ax = axes(use_frame);
    plot(temp)

    title(ax,"Strength of Full Interaction")
    xlabel(ax,"Time (scaled)")
    ylabel(ax,'L^1 of Interaction Matrix')

    saveas(figure(8),'look_at_full_ker.png')

end
