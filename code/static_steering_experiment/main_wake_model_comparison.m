%% Post process the static wake steering experiment results
close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Michael F. Howland, mhowland@mit.edu
% This code is a companion to Howland et al. "Collective wind farm
% operation based on a predictive model increases utility-scale energy 
% production" 
% This script will produce a similar output to Figure 4 in the manuscript,
% although we note that the flow control model predictions here will not be 
% identical to the predictions in the paper because the coefficients of 
% lift, drag, thrust, and power are not available publicly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add paths
addpath(strcat('utility_functions\figplots'))
addpath(strcat('utility_functions'))
addpath(strcat('flow_control_model'))
rng(4096);
goldenrod = [0.85, 0.65, 0.13];

%% Settings
plotNCut = 25; % minimum number of data points to plot a wind direction bin

% Bootstrapping
dataBoot = {};
dataBoot.conf = 0.95; % percent confidence interval
dataBoot.Boots = 1000;

% Wind farm data
load data/turbine_locations % wind farm geometry
load data/lidar_data % LiDAR measurements
load data/turbine_1_yaw % Turbine 1 yaw misalignment values
load data/static_yaw_raw_data % Static yaw misalignment raw data
data=lidar_data;

% LiDAR data
heights = [43,55,67,80,91,104,117,128,141,153,165,200]; 
height_start = 5;
fn = fieldnames(data); 
uv = zeros(length(heights),length(data.(fn{1})));
alpha = uv; w = uv; avail = uv;
for i=1:length(heights);
    uv(i,:) = data.(fn{i+height_start}).('WndSpd');
    alpha(i,:) = data.(fn{i+height_start}).('WndDir');
    avail(i,:) = data.(fn{i+height_start}).('Available');
    w(i,:) = data.(fn{i+height_start}).('Zwind');
end
availVec = any(avail==0,1);

% Load data
TOI = [1, 2, 3, 4];
turbines = {'01','02','03','04'};
rotateTurb = 1; X = X - X(rotateTurb,:); X = X(TOI,:);
% Wind condition bins
directions = [0, 2.5, 5, 7.5, 10, 12.5, 15, 345, 347.5, 350, 352.5, 355, 357.5];
directions_plot = directions; directions_plot(directions_plot>180)=directions_plot(directions_plot>180)-360;
offset_step = 5; offset_bins = -25:offset_step:25;
wind_speed_center = 7; % m/s
% Load wake model data
load data/wakeModelData

% Angle of alignment
alignment_angle = atand( (X(1,2)-X(3,2)) / (X(1,1)-X(3,1)) );
alignment_angle = 360 - (90+alignment_angle);

% Get the data of relevance for each yaw misalignment offset value
for i=1:length(offset_bins);
    yaw{i} = squeeze(output.yaw_store(:,1,i,:));
    Po{i} = squeeze(output.P_opti_store(:,1,i,:));
    Pb{i} = squeeze(output.P_baseline_store(:,1,i,:));
    Pf{i} = squeeze(output.P_fit_store(:,1,i,:));
    Pd{i} = squeeze(output.P_data_store(:,1,i,:));
    Pds{i} = squeeze(output.Prs_data_store(:,1,i,:));
    kw{i} = squeeze(output.kw(:,1,i,:));
    s{i} = squeeze(output.sigma_0(:,1,i,:));
    n{i} = squeeze(output.n(1,i,:));
    Pds_farm{i} = squeeze(output.Pstd_farm_store(1,i,:));
    P_CI_farm{i} = squeeze(output.Pr_farm_CI(:,1,i,:));
    P_CI_all{i} = squeeze(output.P_CI_all_store(:,:,1,i,:));
end

% Initialize wake model
[ turbine, atm, params ] = initalize_inputs(X);
% Wake model setup
params.superposition = 'mod'; 
params.secondary = true;
params_opti.var_p=0.01; params_opti.var_k=5E-5; params_opti.var_sig=5E-5; 
params.Prs = sqrt(params_opti.var_p)*ones(turbine.Nt,1);
params.ucMaxIt = 10; params.Ny = 100; params.epsUc = 10^(-3);
params.yaw_init = zeros(turbine.Nt,1);

%% Power-yaw model inputs
% Below, input the turbine specific airfoil parameters
% The airfoil parameters for the turbine of interest in this experiment are
% not available publicly
airfoil_params = {};
% airfoil_params.R_twist = out.R; airfoil_params.twist = out.twist*pi/180; airfoil_params.c = out.c;
% airfoil_params.interp_setting = out.interp_setting; airfoil_params.airfoil = out.airfoil;
% airfoil_params.interp_pairs = out.interp_pairs; airfoil_params.thickness = out.thickness;
% airfoil_params.tilt = INSERT; airfoil_params.pitch_struct = INSERT;
% airfoil_params.interp_thick = out.interp_thick;
% airfoil_params.tsr = INSERT;
% airfoil_params.zh = INSERT;
% Model settings
params.semi_empirical = false;
cosine_model = [true, true, true, true];
params.local_speed = true;
% Compare to cosine models of the past
params.powerExp=2;

% Initialize loop
Poffset_mat = zeros(turbine.Nt, length(directions), length(offset_bins));
Poffset_std_mat = Poffset_mat; Poffset_store = {};
P_no2nd_offset_mat = zeros(turbine.Nt, length(directions), length(offset_bins));
% Cosine model
Poffset_mat_cos = zeros(turbine.Nt, length(directions), length(offset_bins));
Poffset_std_mat_cos = Poffset_mat; Poffset_store_cos = {};
P_no2nd_offset_mat_cos = zeros(turbine.Nt, length(directions), length(offset_bins));
% Loop over all wind directions
for i=1:length(directions);
    [ XR, indSorted, unsort ] = rotate( X, directions(i), rotateTurb );
    params.kw = kw{offset_bins==0}(indSorted,i); params.sigma_0 = s{offset_bins==0}(indSorted,i);
    turbine.turbCenter = XR; 
    atm.wind_speed = wind_speed_center*ones(turbine.Nt,1); % input real wind speed
    params.cosine_model = cosine_model(indSorted);
    for j = 1:length(offset_bins);
        % Loop over all of the individual 1-min SCADA points and compute
        % the wake model prediction
        Pv = zeros(turbine.Nt, length(output.inds{1,j,i}));
        Pv_cos = zeros(turbine.Nt, length(output.inds{1,j,i}));
        for k=1:length(output.inds{1,j,i});
            local_ind = output.inds{1,j,i}(k);
            local_ind_wind = local_ind;
            % Lidar
            ABL_data = {};
            ABL_data.uv = uv(:, local_ind_wind); 
            ABL_data.alpha = alpha(:, local_ind_wind);
            ABL_data.heights = heights;
            ABL_data.speed_ratio = 0; % INSERT
            if any(isnan(uv(:,output.inds{1,j,i}(k)))==true); 
                % Check for missing LiDAR data and throw away
                Pv(:,k) = NaN;
                Pv_cos(:,k) = NaN;
            else
                % Turbine
                turbine_data = {};
                turbine_data.zhub = 104;
                % Turbine specific details for power-yaw model
                % Yawed turbine
                turbine_data.pitch_yawed = 0; % INSERT
                turbine_data.gamma_set = 0; % INSERT
                turbine_data.lambda_yawed = 0; % INSERT
                turbine_data.WFC_strategy = 0; % INSERT
                % Baseline turbine
                turbine_data.pitch_base = 0; % INSERT
                turbine_data.lambda_base = 0; % INSERT
                turbine_data.gamma_base = 0; % INSERT
                % Store in dictionaries
                turbine.turbine_data = turbine_data;
                turbine.airfoil_params = airfoil_params;
                atm.ABL_data = ABL_data;
                % Run forward model
                turbine.yaw = zeros(turbine.Nt,1);
                turbine.yaw(1) = turbine_1_yaw(local_ind)*pi/180;
                turbine.yaw = turbine.yaw(indSorted);
                [ P, cache, ~ ] = lifting_line_forward_dynamic(turbine, atm, params);
                P = P(unsort); P = P/P(2);
                Pv(:,k) = P;
                % Cosine model comparison
                % Parameters for a cosine model forward pass
                params_cos = params;
                params_cos.powerExp=3;
                params_cos.semi_empirical = false;
                params_cos.local_speed = false;
                params_cos.cosine_model = [true, true, true, true];
                [ Pcos, cache, ~ ] = lifting_line_forward_dynamic(turbine, atm, params_cos);
                Pcos = Pcos(unsort); Pcos = Pcos/Pcos(2);
                Pv_cos(:,k) = Pcos;
            end
        end
        Poffset_mat(:,i,j) = nanmean(Pv,2);
        Poffset_std_mat(:,i,j) = nanstd(Pv,[],2);
        Poffset_store{i,j} = Pv;
        % Cosine model comparison
        Poffset_mat_cos(:,i,j) = nanmean(Pv_cos,2);
        Poffset_std_mat_cos(:,i,j) = nanstd(Pv_cos,[],2);
        Poffset_store_cos{i,j} = Pv_cos;
    end
end

%% Compare power production data to the wake model
for turbines_of_interest = 1; % Option 1: Array power sum. Other options: 1:5;
    switch turbines_of_interest
        case 1
            toi=[1,3,4];
            yb=[0.4,1.6];
            subfolder='turbines_1_3_4/';
        case 2
            toi=[1,3];
            yb=[0.4,1.6];
        case 3
            toi=[1];
            yb=[0.6,1.4];
        case 4
            toi=[3];
            yb=[0,3];
        case 5
            toi=[4];
            yb=[0.4,1.6];
    end
    colors = linspecer(3);
    for dirval = directions(1); % Options: directions
        aind = find(directions==dirval);
        for oind=1:length(offset_bins);
            ptotal(oind) = sum(Pd{oind}(toi,aind));
            pstotal(oind) = sqrt(sum(Pds{oind}(toi,aind).^2));
            ptotal_model(oind) = sum(Poffset_mat(toi,aind,oind));
            ptotal_model_cos(oind) = sum(Poffset_mat_cos(toi,aind,oind));
            pstd_model(oind) = std(sum(Poffset_store{aind,oind}(toi,:),1),[],2);
            ntotal(oind) = n{oind}(aind);
            Pstd_farm(oind) = Pds_farm{oind}(aind);
            Pci_farm(:,oind) = P_CI_farm{oind}(:,aind);
            Pci_all(:,oind) = sqrt(sum(P_CI_all{oind}(:,toi,aind).^2, 2));
            % Bootstrap
            pmodel_sum = sum(Poffset_store{aind,oind}(toi,:),1);
            dataBoot.p = pmodel_sum;
            [stats] = bootstrap_model(dataBoot);
            pCI_model(oind) = stats.P_CI(1);
            % Cut output as per plotNCut value
            if ntotal(oind)<plotNCut;
                ptotal(oind)=NaN; pstotal(oind)=NaN; 
            end
        end
        figure(4); clf; makePlot(4); 
        plot([0,0],yb,'k--', 'LineWidth', 0.75)
        plot([-30,30],[1,1],'k--', 'LineWidth', 0.75);

        % CI from bootstrapping
        p1=shadedErrorBar(offset_bins, ptotal_model/ptotal_model(6), pCI_model/ptotal_model(6), ...
                  'lineProps',{'-','Color', colors(2,:), 'MarkerFaceColor', colors(2,:),'LineWidth',2.75}); 
        p1=p1.mainLine;

        % Cosine model
        green = [0.05, 0.5, 0.06];
        p2=plot(offset_bins, ptotal_model_cos/ptotal_model_cos(6), ...
            '-', 'LineWidth', 1.0, 'color', green);

        if turbines_of_interest==1;
            % CI from bootstrapping
            p3=shadedErrorBar(offset_bins, ptotal/ptotal(6), ...
                Pci_farm/ptotal(6), ...
                   'lineProps',{'d-','Color',colors(1,:),'MarkerFaceColor',colors(1,:),'LineWidth',1.75});
            p3=p3.mainLine;
        else
            shadedErrorBar(offset_bins, ptotal/ptotal(6), Pci_all/ptotal(6), ...
                   'lineProps',{'d-','Color',colors(1,:),'MarkerFaceColor',colors(1,:),'LineWidth',1.75}); 
        end
        xlabel('$\gamma_1$'); ylabel('$\sum^{N_t}_i P_i / \sum^{N_t}_i P_{i}^{\gamma_0}$');
        ylim(yb);
        title(strcat('$\alpha=', num2str(directions_plot(aind)), '^\circ$'))
        
        % Plot model-optimal yaw
        if turbines_of_interest==1 || turbines_of_interest==2;
            [~,ind] = max(ptotal_model); yaw_max = offset_bins(ind);
            p4=plot([yaw_max,yaw_max], yb, 'Color', goldenrod, 'LineWidth', 2);
        end

        % Legend
        if turbines_of_interest==1;
            l = legend([p1,p2,p3,p4],'$\mathrm{Model}$', '$\mathrm{Model},~\hat{P}(\gamma)\sim\cos^3(\gamma)$', ...
                '$\mathrm{Data}$', '$\gamma^*$', 'location', 'southeast');
            set(l,'fontsize',9)
        end
    end
end