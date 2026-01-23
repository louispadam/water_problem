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
    n_omega = parameters.N_x;
    n_theta = parameters.N_y;

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
    data(1,:,:) = ic;
    
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

    Btheta = parameters.B(theta);
    sg_omega = exp(1i*dt*freq_omega.*Btheta');

    ef = angle_to_vec(parameters.ef_angle);
    c = ff_influence(omega);

    for step = 2:steps

        % take Fourier in x
        X_n0 = fft(U,[],2);

        % apply first semigroup
        X_n12 = X_n0 .* sg_omega;

        % invert Fourier in x and compute second semigroup
        U = real(ifft(X_n12,[],2));

        Aomega = parameters.A(omega,ef,U,c);
        sg_theta = exp(1i*dt*Aomega.*freq_theta');

        % apply Fourier in y
        Y_n12 = fft(U,[],1);

        % apply second semigroup
        Y_n1 = Y_n12 .* sg_theta;

        % invert Fourier in y
        U = ifft(Y_n1,[],1);

        tt=tt+dt;

        % Store result at previously calculated frequency
        if step == here
            time(k) = tt;
            data(k,:,:) = real(U);
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

end
