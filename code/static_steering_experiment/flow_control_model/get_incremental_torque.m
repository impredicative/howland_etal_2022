function [dT] = get_incremental_torque(params, u, phi, pitch, r, rho)
% Compute the incremental torque

% Initialize
dT = zeros(length(r), length(phi(1,:)));

% Loop over radial sections
for k=1:length(r);
    % Extract tables
    s(1) = params.interp_pairs(params.interp_setting(k)+1,1);
    s(2) = params.interp_pairs(params.interp_setting(k)+1,2);
    dthick = params.interp_thick(s(2))-params.interp_thick(s(1));
    r1 = (params.interp_thick(s(2))-params.thickness(k))/dthick;
    r2 = (params.thickness(k)-params.interp_thick(s(1)))/dthick;
    aoa_1 = params.airfoil{s(1)}(:,1);
    aoa_2 = params.airfoil{s(2)}(:,1);
    Cl_1 = params.airfoil{s(1)}(:,2);
    Cl_2 = params.airfoil{s(2)}(:,2);
    Cd_1 = params.airfoil{s(1)}(:,3);
    Cd_2 = params.airfoil{s(2)}(:,3);
    % Interp 2
    Cl = r1*Cl_1+r2*Cl_2; Cd = r1*Cd_1+r2*Cd_2;
    Cli = interp1(aoa_1,Cl,phi(k,:)-params.twist(k)-pitch);
    Cdi = interp1(aoa_1,Cd,phi(k,:)-params.twist(k)-pitch);
    % Incremental torque
    dT(k,:) = 0.5*rho*params.c(k).*r(k).*u(k,:).^2 .* ...
        (Cli.*sin(phi(k,:)) - Cdi.*cos(phi(k,:)));
end

end

