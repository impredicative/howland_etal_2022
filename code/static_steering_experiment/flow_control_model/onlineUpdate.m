function [yaws, P_opti, P_baseline, params, error] = onlineUpdate(turbine, atm, params, params_opti)
% Wake steering yaw optimization

% Observe the wind conditions
P_kp1 = params.P; %mean power normalized by first turbine

% State estimation
error = zeros(params_opti.stateEstimationEpochs,1);
psi_kp1 = [repmat(params.kw,[1,params_opti.Ne]); ...
    repmat(params.sigma_0,[1,params_opti.Ne])];
bestParams = params; lowestError = 10; bestStep = 1;
for t=1:params_opti.stateEstimationEpochs;
    [ kw, sigma_0, psi_kp1, ~ ] = EnKF_update( psi_kp1, params_opti.var_k, ...
        params_opti.var_sig, params_opti.var_p, turbine, atm, params, P_kp1 );
    % Metrics
    params.kw = kw; params.sigma_0 = sigma_0;
    [Phat,~] = lifting_line_forward_dynamic(turbine, atm, params);
    % Native power
    %Phat = Phat/10^6;
    % normalized
    Phat = Phat / Phat(1);
    error(t) = mean(abs(Phat-P_kp1)/Phat(1));
    if error(t)<lowestError;
        lowestError=error(t); bestStep = t;
        bestParams=params;
    end
end

% Generate optimal yaw angles
params = bestParams; 
if isfield(params_opti,'uncertain')==0;
    [yaws,P_opti,P_baseline] = yawOptimize(turbine, atm, params, params_opti);
else
    if params_opti.uncertain==1;
        [yaws,P_opti,P_baseline] = yawOptimize_uncertain(turbine, atm, params, params_opti);
    else
        [yaws,P_opti,P_baseline] = yawOptimize(turbine, atm, params, params_opti);
    end
end


end

