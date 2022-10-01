function [yaws,P_opti,P_baseline] = yawOptimize_uncertain(turbine, atm, params, params_opti)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model
% Store original values
X = turbine.turbCenter; turbine.rotateTurb = 1;
[ ~, ~, unsortOrig ] = rotate( X, 270, turbine.rotateTurb );
ws = atm.wind_speed;
paramsFull = params;

% Weibull
if params_opti.weibull==1;
   wsin = params_opti.ws_bins+ws(1);
    params_opti.rho_ws = (params_opti.B/params_opti.A)*...
        (wsin/params_opti.A).^(params_opti.B-1) .* exp(-(wsin/params_opti.A).^params_opti.B);
    params_opti.rho_ws = params_opti.rho_ws / sum(params_opti.rho_ws); 
end
% Baseline
Ptotal = zeros(turbine.Nt,1); yaw_update = turbine.yaw;
for l=1:length(params_opti.p_bins);
    % Get params from one of the power bins
    params = paramsFull{l};
    kw = params.kw(unsortOrig); sigma_0 = params.sigma_0(unsortOrig);
    for j=1:length(params_opti.dir_bins);
        % Rotate turbines for specific wind direction
        [ XR, indSorted, unsort ] = rotate( X, 270+params_opti.dir_bins(j), turbine.rotateTurb );
        turbine.turbCenter = XR; 
        params.kw = kw(indSorted);
        params.sigma_0 = sigma_0(indSorted);
        for i=1:length(params_opti.ws_bins);
            atm.wind_speed = ws + params_opti.ws_bins(i);       
            for mi=1:length(params_opti.gamma_bins);
                turbine.yaw = yaw_update(indSorted);
                % Add offset uncertainty to leading turbine
                turbine.yaw(1) = turbine.yaw(1) + params_opti.gamma_bins(mi) * pi/180;
                % Model
                [ P, cache, ~ ] = lifting_line_forward_dynamic(turbine, atm, params);
                Ptotal = Ptotal + P(unsort) * params_opti.rho_ws(i) * params_opti.rho_dir(j) * ...
                    params.prob(l) * params_opti.rho_gamma(mi);
            end
        end
    end
end
Ptot_baseline = sum(Ptotal); P_baseline = Ptotal;

% Initialize with nonzero yaw?
yaw_update = yaw_update + params.yaw_init;

% eps also determines the termination condition
k=1; check = false; 
m=zeros(turbine.Nt,1); v=zeros(turbine.Nt,1);
bestPower = 0; bestYaw = zeros(turbine.Nt,1); 
bestPowerTurb = zeros(turbine.Nt,1);
Ptot = zeros(params_opti.epochsYaw,1); 
P_time = zeros(params_opti.epochsYaw, turbine.Nt);
yawTime = zeros(params_opti.epochsYaw, turbine.Nt);
P_store = {}; grad_store = {};
while k < params_opti.epochsYaw && check == false;

    % Forward prop
    Ptotal = zeros(turbine.Nt,1); grads_total = Ptotal;
    for l=1:length(params_opti.p_bins);
        % Get params from one of the power bins
        params = paramsFull{l};
        kw = params.kw(unsortOrig); sigma_0 = params.sigma_0(unsortOrig);
        for j=1:length(params_opti.dir_bins);
            % Rotate turbines for specific wind direction
            [ XR, indSorted, unsort ] = rotate( X, 270+params_opti.dir_bins(j), turbine.rotateTurb );
            turbine.turbCenter = XR; 
            params.kw = kw(indSorted);
            params.sigma_0 = sigma_0(indSorted);
            for i=1:length(params_opti.ws_bins);
                atm.wind_speed = ws + params_opti.ws_bins(i);       
                for mi=1:length(params_opti.gamma_bins);
                    turbine.yaw = yaw_update(indSorted);
                    % Add offset uncertainty to leading turbine
                    turbine.yaw(1) = turbine.yaw(1) + params_opti.gamma_bins(mi) * pi/180;
                    % Run model
                    [ P, cache, ~ ] = lifting_line_forward_dynamic(turbine, atm, params);
                    Ptotal = Ptotal + P(unsort) * params_opti.rho_ws(i) * params_opti.rho_dir(j) * ...
                        params.prob(l) * params_opti.rho_gamma(mi);
                    %P_store{i,j,l,m} = P;
                    % Geometry check
                    if isfield(cache, 'Cp')==0;
                        grads.dp_dgamma = 0*grads_total;
                    else
                        [ grads ] = lifting_line_backward_dynamic(turbine, atm, params, cache);
                    end
                    grads_total = grads_total + grads.dp_dgamma(unsort) * ...
                        params_opti.rho_ws(i) * params_opti.rho_dir(j) * ...
                        params.prob(l) * params_opti.rho_gamma(mi);
                    %grads_store{i,j,l,m} = grads.dp_dgamma;
                end
            end
        end
    end
    yawTime(k+1,:) = turbine.yaw;
    Ptot(k) = sum(Ptotal);
    P_time(k,:) = Ptotal;
    yawOld = turbine.yaw;
    % Update best opti result
    if Ptot(k) > bestPower;
       [ ~, indSorted, ~ ] = rotate( X, 270, turbine.rotateTurb );
       bestPower = Ptot(k);
       bestYaw = yaw_update(indSorted);
       bestPowerTurb = Ptotal(indSorted);
    end

    % Gradient ascent update yaw
   if strcmp(params_opti.optimizer, 'vanilla') == true;
        % Vanilla
        turbine.yaw = turbine.yaw + params_opti.learning_rate_yaw*grads.dp_dgamma;
    elseif strcmp(params_opti.optimizer, 'adam') == true;
        % Adam
        m = params_opti.beta1*m + (1-params_opti.beta1)*grads_total;
        v = params_opti.beta2*v + (1-params_opti.beta2)*(grads_total.^2);
        yaw_update = yaw_update + params_opti.learning_rate_yaw .* m ./ (sqrt(v) + params.eps);
    end
    if mod(k, params_opti.printInteger) == 0;
       Ptot(k)/Ptot_baseline 
       %yaw_update(1) * 180/pi
       k
    end
    k = k+1;
    if k > 10;
        if k==20;
            params_opti.learning_rate_yaw = params_opti.learning_rate_yaw / 10;
            m = m / 10; v = v / 10;
        end
        if abs(Ptot(k-1)-Ptot(k-2))/abs(Ptot(k-1)) < 10^-7;%params.eps;
            check = true;
        end
        if isnan(Ptot(k-1))==1 || Ptot_baseline==0;
           check=true; 
        end
    end
end
        
% Final results
yaws = bestYaw * 180/pi;
P_opti = bestPowerTurb;

end

