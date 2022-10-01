function [u] = get_u(params, r, psi, z, alpha_v, uv, w_cond)
%Calculation of u at a fixed (r, theta) location
% Unpack parameters
zh = params.zh;
gamma = params.gamma; 
tsr = params.tsr; R = params.R; a=params.a; ap=1-a;
% Compute 
theta_m = gamma + interp1(z,alpha_v,zh+r.*cos(psi));
theta_m_a = theta_m .* cos(psi);
theta_m_b = theta_m .* sin(psi);
uhub = uv(z==104); 
u = sqrt( ( interp1(z,uv,zh+r.*cos(psi)) .* cos(theta_m_b) .* ...
    cos(theta_m_a) * ap).^2 + ...
    ( uhub.*tsr*r*w_cond/R - ...
    interp1(z,uv,zh+r.*cos(psi)) .* cos(theta_m_b) .* sin(theta_m_a) * ap ).^2 );
end

