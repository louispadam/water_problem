function return_data = A(z,curr_state,conv_data,params)
%A Calculation of acceleration, which includes an interaction (so
%convolution) and an electric field. The nomenclature is an artifact of the
%more general derivation of pseudospectral-etd methods.
%
%last updated 03/10/26 by Adam Petrucci

    % retrieve electric field
    e_vec = angle_to_vec(params.ef_angle);

    % integrate omega out of state (acceleration depends only on theta)
    no_omega = squeeze(trapz(curr_state,2));
    atv = angle_to_vec(z);
    no_omega_theta = squeeze(trapz(no_omega .* permute(atv, [3 2 1]),2));

    % compute interaction by convolution against given kernel
    interaction_1 = conv(squeeze(conv_data(:,1)),squeeze(no_omega_theta(:,1)),"same") + ...
                    conv(squeeze(conv_data(:,2)),squeeze(no_omega_theta(:,2)),"same");
    interaction_2 = conv(squeeze(conv_data(:,3)),squeeze(no_omega_theta(:,1)),"same") + ...
                    conv(squeeze(conv_data(:,4)),squeeze(no_omega_theta(:,2)),"same");
    interaction_total = [interaction_1, interaction_2];
    
    % incorporate electic field
    net = params.ker_pow * interaction_total - ...
            params.ef_pow * e_vec';

    return_data = net * angle_to_perp(z);

end