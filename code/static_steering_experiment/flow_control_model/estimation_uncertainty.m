function [bestParams, error] = estimation_uncertainty(turbine, atm, params, params_opti)
% Predicts the power production given a yaw misalignment strategy and
% baseline power production

% Initialize the uncertain power production
bestParams = {};
stats = {}; stats.P = mean(params.P); stats.Prs = params.Prs;
paramsOrig = params; paramsOrig.Phat = zeros(turbine.Nt,1);

for i = 1:length(params_opti.p_bins);
    % Initialize
    params = paramsOrig;
    params.seed = params.seed + params_opti.stateEstimationEpochs * params_opti.p_bins(i);
    params.P = params.P + params_opti.p_bins(i)*params.Prs;
    params.P(params.P>1) = 1;
    % Observe the wind conditions
    P_kp1 = params.P; %mean power normalized by first turbine
    params.prob(i) = params_opti.prob(i);
    % State estimation
    error = zeros(params_opti.stateEstimationEpochs,1);
    psi_kp1 = [repmat(params.kw,[1,params_opti.Ne]); ...
        repmat(params.sigma_0,[1,params_opti.Ne])];
    bestParams{i} = params; lowestError = 10; bestStep = 1;
    for t=1:params_opti.stateEstimationEpochs;
        params.seed = params.seed+1;
        [ kw, sigma_0, psi_kp1, ~ ] = EnKF_update( psi_kp1, params_opti.var_k, ...
            params_opti.var_sig, params_opti.var_p, turbine, atm, params, P_kp1 );
        % Metrics
        params.kw = kw; params.sigma_0 = sigma_0;
        [Phat,~] = lifting_line_forward_dynamic(turbine, atm, params);
        % Native power
        %Phat = Phat/10^6;
        % normalized
        Phat = Phat / Phat(turbine.normalizer);
        params.Phat = Phat;
        error(t) = mean(abs(Phat-P_kp1)/Phat(turbine.normalizer));
        if error(t)<lowestError;
            lowestError=error(t); bestStep = t;
            bestParams{i}=params;
        end
    end
end

end

