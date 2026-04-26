function return_data = A(z,curr_state,conv_data,mass,params)
%A Calculation of acceleration, which includes an interaction (convolution
%and possibly point mass) and an electric field. The nomenclature is an
%artifact of the more general derivation of pseudospectral-etd methods.
%
%last updated 03/16/26 by Adam Petrucci

    % retrieve electric field
    e_vec = angle_to_vec(params.ef_angle);

    % integrate omega out of state (acceleration depends only on theta)
    no_omega = squeeze(trapz(curr_state,2));
    atv = angle_to_vec(z);
    no_omega_theta = squeeze(trapz(no_omega .* permute(atv, [3 2 1]),2));

    %no_omega_theta_per = [no_omega_theta, no_omega_theta, no_omega_theta];
    %interaction_1 = conv(squeeze(conv_data(:,1)),squeeze(no_omega_theta_per(:,1)),"same") + ...
    %                conv(squeeze(conv_data(:,2)),squeeze(no_omega_theta_per(:,2)),"same");
    %interaction_2 = conv(squeeze(conv_data(:,3)),squeeze(no_omega_theta_per(:,1)),"same") + ...
    %                conv(squeeze(conv_data(:,4)),squeeze(no_omega_theta_per(:,2)),"same");
    %interaction_total = [interaction_1, interaction_2] + (mass * no_omega_theta')';

    interaction_1f = fft(squeeze(conv_data(:,1))) .* fft(ifftshift(squeeze(no_omega_theta(:,1)))) + ...
                     fft(squeeze(conv_data(:,2))) .* fft(ifftshift(squeeze(no_omega_theta(:,2))));
    interaction_2f = fft(squeeze(conv_data(:,3))) .* fft(ifftshift(squeeze(no_omega_theta(:,1)))) + ...
                     fft(squeeze(conv_data(:,4))) .* fft(ifftshift(squeeze(no_omega_theta(:,2))));
    interaction_1f = ifft(interaction_1f);
    interaction_2f = ifft(interaction_2f);
    %if mod(length(interaction_1f),2) == 0
    %    interaction_1f(1) = (interaction_1f(1) + interaction_1f(end))/2;
    %    interaction_1f(2:end) = (interaction_1f(2:end) + interaction_1f(1:end-1))/2;
    %    interaction_2f(1) = (interaction_2f(1) + interaction_2f(end))/2;
    %    interaction_2f(2:end) = (interaction_2f(2:end) + interaction_2f(1:end-1))/2;
    %end
    interaction_total = [interaction_1f, interaction_2f] + (mass * no_omega_theta')';
    
    % incorporate electric field
    net = params.ker_pow * interaction_total - params.ef_pow * e_vec';

    return_data = net * angle_to_perp(z);

end