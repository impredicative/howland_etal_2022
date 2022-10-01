function [ powerOut, CpOut, CtOut, dCt_dyaw, eta ] = lookup_tables(ws, A, yaw, atm, two);
%Experimental Ct and Cp look up tables
% Parameters
rho = atm.rho;

% 5-25 m/s
wind_speed = linspace(3, 20, 20-2);
wind_speed_ct = wind_speed;

% Power curve
poly_1 = 1.0e+04 * [0.303351137777127  -4.218819053193799];
poly_2 = 1.0e+06 * [0.000616236162362   1.459763837638376];
p_model = poly_1(1)*wind_speed.^3 + poly_1(2);
ind=find(wind_speed>=9 & wind_speed<11);
p_model(ind) = poly_2(1)*wind_speed(ind).^3 + poly_2(2);
ind=find(wind_speed>=11);
p_model(ind)=2100*1000;
power = p_model;

% Thrust curve
Ct_model = zeros(length(wind_speed),1);
Ct_model(wind_speed<8) = 0.83;
ind = find(wind_speed>=8);
poly = [-0.000714160839161   0.037802947052947  -0.675690934065938   4.160945804195820];
Ct_model(ind) = polyval(poly,wind_speed(ind));
Ct = Ct_model;

% Out values
if ws >= wind_speed(1) && ws <= wind_speed(end);
    powerOut = interp1(wind_speed, power, ws);
    CpOut = powerOut / (0.5*rho*A*ws^3);
    CtOut = interp1(wind_speed_ct, Ct, ws); 
    ap = 0.5*(1-sqrt(1-CtOut));
    eta = CpOut / (4*ap*(1-ap)^2);
elseif ws > wind_speed(end);
    powerOut = power(end);
    CpOut = powerOut / (0.5*rho*A*ws^3);
    CtOut = Ct(end);
    ap = 0.5*(1-sqrt(1-CtOut));
    eta = CpOut / (4*ap*(1-ap)^2);
elseif ws < wind_speed(1); % below or above the cut-in speed
    powerOut = 0;
    CpOut = 0; 
    CtOut = 0; 
    eta = 0;
end

dCt_dyaw = 0;

end

