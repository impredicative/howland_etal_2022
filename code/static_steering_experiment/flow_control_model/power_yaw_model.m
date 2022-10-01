function [Pr] = power_yaw_model(ABL_data, turbine_data, params, semi_empirical)
% Compute the power production for a wind turbine based on the incident
% wind field and the yaw misalignment of the turbine
% See model details in Howland et al. JRSE (2020)

% Initialize model
nt = 50; % number of points in theta direction
rho = 1.225;
r = params.R_twist;
theta = linspace(0,2*pi,nt);
params.R = params.R_twist(end);
rmat = repmat(r,[1,length(theta)]);
thetamat = repmat(theta,[length(r),1]);

% Get LiDAR conditions
alpha_l = ABL_data.alpha - ABL_data.alpha(ABL_data.heights==turbine_data.zhub);
alpha_l(alpha_l<-180)=alpha_l(alpha_l<-180)+360;
alpha_l(alpha_l>180)=alpha_l(alpha_l>180)-360;
u_l = ABL_data.uv;

% Compute the angular velocity of the yawed turbine
params.gamma = turbine_data.gamma_yawed; 
if semi_empirical==true; % Uses the measured angular velocity 
    w_norm = turbine_data.lambda_yawed / params.tsr;
    params.wb = turbine_data.lambda_base / params.tsr;
    Or_model = w_norm / params.wb;
else % Estimates the angular velocity by solving an optimization problem
    params.heights=ABL_data.heights; params.alpha_l=alpha_l; params.a = 1/3;
    params.u_l=u_l; params.theta = theta; params.r = r;
    params.pitch = turbine_data.pitch_yawed; 
    params.yaw_base = turbine_data.gamma_base;
    params.pitch_base = turbine_data.pitch_base;
    params.tsr = turbine_data.lambda_base;
    params.wb = 1;
    gamma_corrected = turbine_data.gamma_set;
    w_norm = 1;
    Or_model = w_norm / params.wb;
    end
end

%% Run forward pass of model to compute power-yaw relationship
params.a = 1/3;
u = get_u(params, rmat, thetamat, ABL_data.heights, alpha_l*pi/180, u_l, w_norm);
phi = get_phi(params, rmat, thetamat, ABL_data.heights, alpha_l*pi/180, u_l, w_norm);
dT = get_incremental_torque(params, u, phi, turbine_data.pitch_yawed, r, rho);

%% Aligned turbine power
params.a = 1/3;
params.gamma = turbine_data.gamma_base; 
u = get_u(params, rmat, thetamat, ABL_data.heights, alpha_l*pi/180, u_l, params.wb);
phi = get_phi(params, rmat, thetamat, ABL_data.heights, alpha_l*pi/180, u_l, params.wb);
dTa = get_incremental_torque(params, u, phi, turbine_data.pitch_base, r, rho);

% Omega ratio
Or = Or_model;
% Torque ratio
T = trapz(theta, trapz(r, dT, 1), 2) / (trapz(theta, trapz(r, dTa, 1), 2));
% Power ratio
Pr = T*Or; 

end

