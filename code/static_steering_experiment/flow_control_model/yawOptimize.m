function [yaws,P_opti,P_baseline] = yawOptimize(turbine, atm, params, params_opti)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model
% Load model inputs
[ P, cache, ~ ] = lifting_line_forward_dynamic(turbine, atm, params);
Ptot_baseline = sum(P); P_baseline = P;

% Initialize with nonzero yaw?
turbine.yaw = turbine.yaw + params.yaw_init;

% eps also determines the termination condition
k=1; check = false; 
m=zeros(turbine.Nt,1); v=zeros(turbine.Nt,1);
bestPower = 0; bestYaw = zeros(turbine.Nt,1); 
bestPowerTurb = zeros(turbine.Nt,1);
Ptot = zeros(params_opti.epochsYaw,1); 
P_time = zeros(params_opti.epochsYaw, turbine.Nt);
yawTime = zeros(params_opti.epochsYaw, turbine.Nt);
while k < params_opti.epochsYaw && check == false;

    % Forward prop
    [ P, cache, ~ ] = lifting_line_forward_dynamic(turbine, atm, params);
    yawTime(k+1,:) = turbine.yaw;
    Ptot(k) = sum(P);
    P_time(k,:) = P;
    yawOld = turbine.yaw;
    if isfield(cache, 'Cp')==0;
        k=k+1;
        bestPowerTurb = P;
        break
    end
    [ grads ] = lifting_line_backward_dynamic(turbine, atm, params, cache);
    % Update best opti result
    if Ptot(k) > bestPower;
       bestPower = Ptot(k);
       bestYaw = turbine.yaw;
       bestPowerTurb = P;
    end

    % Gradient ascent update yaw
   if strcmp(params_opti.optimizer, 'vanilla') == true;
        % Vanilla
        turbine.yaw = turbine.yaw + params_opti.learning_rate_yaw*grads.dp_dgamma;
    elseif strcmp(params_opti.optimizer, 'adam') == true;
        % Adam
        m = params_opti.beta1*m + (1-params_opti.beta1)*grads.dp_dgamma;
        v = params_opti.beta2*v + (1-params_opti.beta2)*(grads.dp_dgamma.^2);
        turbine.yaw = turbine.yaw + params_opti.learning_rate_yaw .* m ./ (sqrt(v) + params.eps);
    end
    if mod(k, params_opti.printInteger) == 0;
       Ptot(k)/Ptot_baseline 
       k
    end
    k = k+1;
    if k > 10;
%         if abs(Ptot(k-1)-Ptot(k-2))/abs(Ptot(k-1)) < 10^-16;%params.eps;
%             check = true;
%         end
        if isnan(Ptot(k-1))==1 || Ptot_baseline==0;
           check=true; 
        end
    end
end

% Final results
yaws = bestYaw * 180/pi;
P_opti = bestPowerTurb;

end

