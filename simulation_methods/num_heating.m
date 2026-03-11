function return_data = num_heating()
y_heat     = ones([1,parameters.N_y]);
omega_heat = abs(omega)>0.9*pi;
theta_heat = ones([1,parameters.N_theta]);

cut_heat = y_heat' .* omega_heat .* reshape(theta_heat,1,1,[]);

lt = length(t_c);

grain = 100;

track_heat = zeros(1,grain);

for i = 1:grain
    j = round(i*lt/grain);
    track_heat(i) = trapz(y, trapz(omega, trapz(theta, cut_heat .* abs(squeeze(d_c(j,:,:,:))), 3), 2));
    if mod(i,round(grain/10)) == 0
        fprintf('%d  Simulation Progress: %3.0f%%\n',i,100*i/grain);
    end
end
return_data = [t_c(round((1:grain)*lt/grain));track_heat];
end