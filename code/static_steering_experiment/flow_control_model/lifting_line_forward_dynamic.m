function [ P, cache, uc_error ] = lifting_line_forward_dynamic(turbine, atm, params);
%Implement the forward propagation of the lifting line model
% Inputs: domain dictionary
% Outputs: power per turbine (P), cache required for backprop (cache)
% See: 
% Relevant papers:
% 1) Howland, Michael F., Sanjiva K. Lele, and John O. Dabiri. 
% "Wind farm power optimization through wake steering." 
% Proceedings of the National Academy of Sciences 
% 116.29 (2019): 14495-14500.
% 2) Howland, M. F., Ghate, A. S., Lele, S. K., and Dabiri, J. O.: 
% Optimal closed-loop wake steering, Part 1: Conventionally neutral 
% atmospheric boundary layer conditions, 
% Wind Energ. Sci., https://doi.org/10.5194/wes-2020-52, 2020. 
% 3) Howland, Michael F., and John O. Dabiri. "Influence of wake model 
% superposition and secondary steering on model-based wake steering control 
% with SCADA data assimilation." Energies 14, no. 1 (2020): 52.
% 4) Howland, Michael F., et al. "Collective wind farm operation based on a 
% predictive model increases utility-scale energy production." 
% arXiv preprint arXiv:2202.06683 (2022).

%%% Initialize
%%%%%%%%%%%%%%%% Do not change anything between these lines %%%%%%%%%%%%%%%
% Unpack from dictionaries
% Turbine info
Nt = turbine.Nt; turbCenter = turbine.turbCenter;
isFront = turbine.isFront; isBack = turbine.isBack; %graph
%leadingTurbine = turbine.leadingTurbine; 
%followingTurbine = turbine.followingTurbine;
D = turbine.D; D = ones(Nt,1) * D;
A = D.^2 * pi/4;
D = D ./ D; R = D/2; % by convention D=1 
yaw = turbine.yaw; % * pi / 180; % convert to radians
% Atmospheric conditions
rho = atm.rho; wind_speed = atm.wind_speed;
% Model parameters
kw = params.kw; sigma_0 = params.sigma_0;
delta_w = R; eps = params.eps;
Nx = params.Nx; % for spatial integration
% Initialize variables in function
cache = {};
uInf = zeros(Nt,1); delta_v0 = zeros(Nt,1); a = zeros(Nt,1); ap = zeros(Nt,1);
delta_u0 = zeros(Nt,1); u_eff = zeros(Nt,1); Cp = zeros(Nt,1);
P = zeros(Nt,1); %yCenter = zeros(Nt,1); 
%sigmaEnd = zeros(Nt,1);
y_c_end = zeros(Nt,1);
duCache = zeros(Nt,1); du_no_a_cache = zeros(Nt,1);
duSquareSum = zeros(Nt,1); gaussianStore = zeros(Nt,1);
dCt_dyaw = zeros(Nt,1);
Uc = ones(Nt,1)*wind_speed(1); Uc_old = Uc; ucErrorStore = 0;
uci = ones(Nt,Nt)*wind_speed(1); ucr = ones(Nt,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forward prop
% Loop over turbines
uc_error = Inf; ucCount = 1;
while uc_error>params.epsUc && ucCount<params.ucMaxIt;
    % Initialization
    v = zeros(Nt,1);
    for j = 1:Nt;
        % Find the leading turbine
        turbine.turbinesInFront{j} = [];
        for i = 1:Nt;
            dx = turbCenter(j, 1) - turbCenter(i, 1);
            if dx > 0;
               dw = 1 + kw(i) * log(1 + exp((dx - 2.0 .* delta_w(i)) / R(i)));
               boundLow = turbCenter(i,2)-dw/2;
               boundHigh = turbCenter(i,2)+dw/2;
               edgeLow = turbCenter(j,2)-D(j)/2;
               edgeHigh = turbCenter(j,2)+D(j)/2;
               if (edgeLow>=boundLow & edgeLow<=boundHigh) || ...
                       (edgeHigh>=boundLow & edgeHigh<=boundHigh);
                  turbine.turbinesInFront{j} = [turbine.turbinesInFront{j}, i]; 
               end
            end
        end
        % Loop over front turbines
        delta_u_face = 0; duSum = zeros(params.Ny,1); us = duSum;
        ylocal = linspace(turbCenter(j, 2)-D(j)/2, turbCenter(j, 2)+D(j)/2, params.Ny)';
        for k = turbine.turbinesInFront{j};
            aPrev = 0.5*(1-sqrt(1-cache.Ct(k)*cos(yaw(k))^2));
            dx = turbCenter(j, 1) - turbCenter(k, 1);
            xpFront = linspace(0, dx, Nx);
            dw = 1 + kw(k) * log(1 + exp((xpFront - 2.0 .* delta_w(k)) / R(k)));
            if strcmp(params.superposition, 'mom') == true || ...
                    strcmp(params.superposition, 'mod') == true;
                du = u_eff(k) * aPrev/(dw(end)^2+params.eps) * (1+erf(xpFront(end)/(delta_w(k)*sqrt(2))));
            else
                du = wind_speed(j) * aPrev/(dw(end)^2+params.eps) * (1+erf(xpFront(end)/(delta_w(k)*sqrt(2))));
            end
            deltaUIndividual{j}(k) = du;
            duSquareSum(j) = du;

            % Deflection
            dx = turbCenter(j, 1) - turbCenter(k, 1);
            % Secondary steering
            if strcmp(params.superposition, 'linear') == true;
                dvi = wind_speed(j) * 0.25 * cache.Ct(k) * cos(yaw(k))^2 * sin(yaw(k));
            elseif strcmp(params.superposition, 'mom') == true || ...
                    strcmp(params.superposition, 'mod') == true;
                dvi = u_eff(k) * 0.25 * cache.Ct(k) * cos(yaw(k))^2 * sin(yaw(k));
            end
            if params.secondary==true;
                delta_v0(j) = dvi + v(k);
            else
                delta_v0(j) = dvi;
            end
            % Local x frame
            xp{k, j} = linspace(0, dx, Nx);
            % Velocity deficit
            dw = 1 + kw(k) * log(1 + exp((xp{k, j} - 2.0 .* delta_w(j)) / R(j)));
            delta_v = (delta_v0(j) ./ (dw.^2+params.eps)) * 0.5 .* (1 + erf(xp{k, j} ./ (delta_w(j)*sqrt(2))));
            delta_v(xp{k, j}<0) = 0;
            dv = (1 ./ (dw(end)^2+params.eps)) * 0.5 .* (1 + erf(dx ./ (delta_w(j)*sqrt(2))));
            dvec = (1 ./ (dw.^2+params.eps)) * 0.5 .* (1 + erf(xp{k, j} ./ (delta_w(j)*sqrt(2))));
            % Secondary steering
            if strcmp(params.superposition, 'linear') == true;
                y_c{k, j} = cumtrapz(xp{k, j}, -delta_v/wind_speed(j));
            elseif strcmp(params.superposition, 'mom') == true || ...
                    strcmp(params.superposition, 'mod') == true;
                y_c{k, j} = cumtrapz(xp{k, j}, -delta_v/u_eff(k));
            end
            %Local y frame
            sigma = sigma_0(k) * dw; 
            sigmaEnd(k, j) = sigma(Nx); % take the last value of sigma which is at the next turbine
            ind = Nx;
            yCenter(k, j) = y_c{k, j}(ind) + turbCenter(k, 2);
            y_c_end(k, j) = y_c{k, j}(ind);
            % Effective velocity at turbine face
            gaussian = erf((turbCenter(j, 2)+D(j)/2-yCenter(k, j))/sqrt(2*sigmaEnd(k, j)^2)) - ...
                erf((turbCenter(j, 2)-D(j)/2-yCenter(k, j))/sqrt(2*sigmaEnd(k, j)^2));
            gaussian = gaussian * (D(k)^2)/(16*sigma_0(k)^2) * sigmaEnd(k, j)*sqrt(2*pi);
            if strcmp(params.superposition, 'sos') == true;
                delta_u_face = delta_u_face + (du * gaussian)^2;
            elseif strcmp(params.superposition, 'mom') == true;
                % compute local convective velocity
                us = du * D(k)^2 / (8*sigma_0(k)^2) * exp(-(ylocal-yCenter(k, j)).^2 / (2*sigmaEnd(k, j)^2));
                uw = u_eff(k) - us;
                uci(k, j) = trapz(ylocal, uw.*us) / trapz(ylocal, us);
                % Update cache stuff
                delta_u_face = delta_u_face + du * gaussian * uci(k, j) / Uc(j);
                duSum = duSum + us * uci(k, j) / Uc(j);
                deltaUIndividual{j}(k) = du * uci(k, j) / Uc(j);
                erfStore(k, j) = (du/u_eff(k)) * gaussian * uci(k, j) / Uc(j);
            elseif strcmp(params.superposition, 'mod') == true;
                uci(k, j) = wind_speed(1);
                % Update cache stuff
                delta_u_face = delta_u_face + du * gaussian * uci(k, j) / Uc(j);
                duSum = duSum + us * uci(k, j) / Uc(j);
                deltaUIndividual{j}(k) = du * uci(k, j) / Uc(j);
                erfStore(k, j) = (du/u_eff(k)) * gaussian * uci(k, j) / Uc(j);
            else
                delta_u_face = delta_u_face + du * gaussian;
            end
            gaussianStore(k, j) = gaussian;
            % Secondary steering lateral velocity at the turbine face
            if params.secondary == true;
                ucr(k, j) = (uci(k, j) / Uc(j)) * dv * gaussian;
                v(j) = v(j) + dvi * ucr(k, j);
            end
        end
        if strcmp(params.superposition, 'sos') == true;
            duSquareSum(j) = delta_u_face;
            delta_u_face = sqrt(delta_u_face);
        elseif strcmp(params.superposition, 'mom') == true;
            % Global convection velocity
            Uc_old(j) = Uc(j); us = duSum;
            uw = wind_speed(j) - us;
            Uc(j) = trapz(ylocal, uw.*us) / trapz(ylocal, us);
            if isnan(Uc(j))==true;
                Uc(j) = wind_speed(j);
            end
        elseif strcmp(params.superposition, 'mod') == true;
            Uc_old(j) = wind_speed(1); Uc(j) = Uc(j);
        end
        delta_u_face_store(j) = delta_u_face;
        % Effective velocity
        u_eff(j) = (1/D(j)) * (D(j)*wind_speed(j) - delta_u_face);
        if u_eff(j) >= 0;
            [ ~, ~, cache.Ct(j), dCt_dyaw(j), cache.eta(j)] = lookup_tables(u_eff(j), A(j), 0, atm, turbine.turbine_type{j});
        else
            u_eff(j) = 0;
            [ ~, ~, cache.Ct(j), dCt_dyaw(j), cache.eta(j)] = lookup_tables(u_eff(j), A(j), 0, atm, turbine.turbine_type{j});
        end
        % Calculate the induction factor at the end
        a(j) = 0.5*(1-sqrt(1-cache.Ct(j)*cos(yaw(j))^2));
        % Induction factor for power
        ap(j) = 0.5*(1-sqrt(1-cache.Ct(j)));
        % Power-yaw model
        % Cosine model or Howland et al. JRSE (2020) model
        % CAUTION: Only use Howland model for freestream turbines
        % since we don't know the velocity profile incident to downwind
        % turbines under waked inflow from a yawed turbine
        if params.cosine_model(j)==false;
            turbine.turbine_data.gamma_yawed = yaw(j);
            Pr = power_yaw_model(atm.ABL_data, turbine.turbine_data, ...
                turbine.airfoil_params, params.semi_empirical);
            Cp(j) = 4*cache.eta(j)*ap(j)*(1-ap(j))^2 * Pr;
        else
            Cp(j) = 4*cache.eta(j)*ap(j)*(1-ap(j))^2 * cos(yaw(j))^params.powerExp;
        end
        % Compute power
        if u_eff(j) < 25;
            P(j) = 0.5 * rho * A(j) * Cp(j) * u_eff(j)^3;  
        else
            P(j) = params.ratedPower; % rated wind speed
        end
    end
    % Error update
    if strcmp(params.superposition, 'mom') == true;
        uc_error = mean(abs(Uc_old-Uc)); ucErrorStore(ucCount) = uc_error;
        ucCount = ucCount+1;
    else
        uc_error=0;
    end
end

% Store some stuff in the cache
if exist('xp')==1;
    cache.Cp = Cp; cache.u_eff = u_eff; cache.a = a; cache.xp = xp;
    cache.yCenter = yCenter; cache.sigmaEnd = sigmaEnd;
    cache.uInf = uInf; cache.ap = ap; cache.y_c = y_c_end;
    cache.deltaUIndividual = deltaUIndividual;
    cache.turbinesInFront = turbine.turbinesInFront;
    cache.delta_u_face_store = delta_u_face_store;
    cache.duSquareSum = duSquareSum;
    cache.gaussianStore = gaussianStore; cache.ucErrorStore = ucErrorStore;
    cache.turbinesInFront = turbine.turbinesInFront;
    if strcmp(params.superposition, 'mom') == true || ...
            strcmp(params.superposition, 'mod') == true;
        cache.erfStore = erfStore; cache.uci = uci; cache.Uc = Uc;
        cache.ucr = ucr; cache.v = v;
    end
end

end

