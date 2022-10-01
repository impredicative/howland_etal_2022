function [ kw, sigma_0, psi_kp1, stats ] = EnKF_update( psi_k, var_k, var_sig, var_p, turbine, atm, params, P_kp1 )
% EnKF update
% Closed-loop wake steering control under uncertainty
% Relevant papers:
% 1) Howland arXiv 2021
% 2) Howland, M. F., Ghate, A. S., Lele, S. K., and Dabiri, J. O.: 
% Optimal closed-loop wake steering, Part 1: Conventionally neutral 
% atmospheric boundary layer conditions, 
% Wind Energ. Sci. Discuss., https://doi.org/10.5194/wes-2020-52, in review, 2020. 
rng(params.seed);

% Define relevant parameters
[NN, Ne] = size(psi_k); Nt = NN/2;

% Intermediate forcast step
chi = [randn(Nt,Ne)*sqrt(var_k); randn(Nt,Ne)*sqrt(var_sig)];
I = eye(2*Nt,2*Nt); psi_kp = zeros(size(psi_k));
psiHat_kp = zeros(Nt,Ne);
psi_kp = psi_k + I*chi;
kwPerturb = psi_kp(1:Nt,:);
sigmaPerturb = psi_kp(Nt+1:2*Nt,:);
% Zero yaw normalization
turbine_zeroYaw = turbine; turbine_zeroYaw.yaw = turbine_zeroYaw.yaw*0;
[Phat_zeroYaw,~] = lifting_line_forward_dynamic(turbine_zeroYaw, atm, params);
for i=1:Ne; % can be parfor
    % Lifting line model
    paramsLocal = params;
    paramsLocal.kw = kwPerturb(:,i);
    paramsLocal.sigma_0 = sigmaPerturb(:,i);
    [Phat,~] = lifting_line_forward_dynamic(turbine, atm, paramsLocal);
    % Native power
    %psiHat_kp(:,i) = Phat/10^6;
    % Normalize by first turbine
    psiHat_kp(:,i) = Phat/Phat_zeroYaw(1);
end


% Perturbation
one_Ne = ones(Ne, Ne) / Ne;
psi_mean = psi_kp * one_Ne;
psiHat_mean = psiHat_kp * one_Ne;
psi_kpp = psi_kp - psi_mean;
psiHat_kpp = psiHat_kp - psiHat_mean;

% Perturbation matrices
E_kp1 = randn(Nt, Ne) * sqrt(var_p);
burger_kp1 = repmat(P_kp1, [1, Ne]) + E_kp1;


% Measurement analysis step
psi_kp1 = psi_kp + psi_kpp * psiHat_kpp' * inv(psiHat_kpp*psiHat_kpp' + ...
    E_kp1 * E_kp1') * (burger_kp1 - psiHat_kp);

% Ouput the final values
psi_kp1_mean = psi_kp1*one_Ne;
kw = psi_kp1_mean(1:Nt,1);
sigma_0 = psi_kp1_mean(Nt+1:2*Nt,1);
error = mean(abs(psiHat_mean(:,1)-P_kp1));
meanOut = psiHat_mean(:,1);
stdOut = std(psiHat_kp,0,2);
stats = {};
stats.meanOut = meanOut; stats.stdOut = stdOut;

end

