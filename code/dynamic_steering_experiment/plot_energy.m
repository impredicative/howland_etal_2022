%% Plot farm energy production

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Michael F. Howland, mhowland@mit.edu
% This code is a companion to Howland et al. "Collective wind farm
% operation based on a predictive model increases utility-scale energy 
% production" 
% This script will produce a similar output to Figure 5 in the manuscript,
% except the energy ratio weighting uses the 'switched' methdology here
% (see Supplementary Note 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

% Add paths
addpath(strcat('utility_functions\figplots'))
addpath(strcat('utility_functions'))

%% Lookup table
load('data/lookup_tables.mat');

% Load the energy gain data
load('data/data_6_8') % Wind speeds between 6 and 8 m/s
load('data/data_0_20') % Wind speeds between 0 and 20 m/s

% Load data count
load('data/N_6_8') % Histogram data 6-8 m/s
load('data/N_0_20') % Histogram data 0-20 m/s

% Wind direction bins
directions = 0:2.5:360;

% Load sectors for energy gain calculations
[E_dir_bounds] = E_gain_sectors(directions);

% Bounds for the plot
xb = [-20, 15]; 
inflection = -5;
xbp = [inflection-20, inflection+20];
yb = [0.25,1.2]; yb2 = [-0.2,0.2];
[dirw,indw] = sort(wrapTo180(directions));

% Include histogram as the subplot
colors = linspecer(5); indc = 4;
f6=figure(6); clf; makePlot(6); hold on
set(f6, 'Position', [50 50 500 375] )
sp1=subplot(3,1,1); hold on; box on;
plot(xbp,[0,0],'k-','LineWidth',1.75);
plot([inflection, inflection],yb2*100,'k-','LineWidth',1.75);
plot([xb(1),xb(1)],yb2*100,'k--','LineWidth',1.75);
plot([xb(2),xb(2)],yb2*100,'k--','LineWidth',1.75);
% Plot
% Full wind speed range plot
shadedErrorBar(dirw, out_0_20.E_gain_y_0(indw,1,2)*100, out_0_20.E_gain_y_0_CI(:,indw,1,2)*100, ...
    'lineProps',{'o--','Color',colors(1,:),'MarkerFaceColor',colors(1,:),'LineWidth',2}); 
% Low wind speed plot
shadedErrorBar(dirw, out_6_8.E_gain_y_0(indw,1,2)*100, out_6_8.E_gain_y_0_CI(:,indw,1,2)*100, ...
    'lineProps',{'o--','Color',colors(indc,:),'MarkerFaceColor',colors(indc,:),'LineWidth',2}); 
% Aesthetics
plot([xb(1),xb(1)],yb2*100,'k--','LineWidth',1.75);
plot([xb(2),xb(2)],yb2*100,'k--','LineWidth',1.75);
xlim([xb(1)-5, xb(2)+5]); ylim([-10,20]); 
ylabel('$\mathrm{Energy \ gain} \ (\%)$');
pos1 = get(sp1, 'Position');
pos1_mod = pos1 + [0.01 0 0 0.05];
set(sp1, 'Position', pos1_mod )
grid on;
set(gca,'Xticklabel',[])

sp2 = subplot(3,1,2); hold on; grid on;
[dirw_bar,indw_bar] = sort(wrapTo180(directions(2:end)));
Nbar_0_20 = squeeze(sum(N_0_20(2:end,:,2,2), 2));
ind = find(dirw_bar>wrapTo360(xb(1)-5) | dirw_bar<wrapTo360(xb(2)+5));
bar(dirw_bar, Nbar_0_20(indw_bar));
% Low wind speed plot
Nbar_6_8 = squeeze(sum(N_6_8(2:end,:,2,2), 2));
bar(dirw_bar, Nbar_6_8(indw_bar), 'FaceColor', colors(indc,:));
plot([xb(1),xb(1)],[-10000,10000],'k--','LineWidth',1.75);
plot([xb(2),xb(2)],[-10000,10000],'k--','LineWidth',1.75);
xlim([xb(1)-5, xb(2)+5]);
ylim([0,1400])
pos2 = get(sp2, 'Position');
pos2_mod = pos2 - [-0.01 -0.03 0 0.05];
set(sp2, 'Position', pos2_mod )
ylabel('$\mathrm{Samples}$');
legend('$0$\textless$u$\textless$20$ m/s', '$6$\textless$u$\textless$8$ m/s', 'location', 'northwest','numcolumns',2);
% Aesthetics
set(sp2,'DefaultAxesFontSize',18);
set(sp2,'DefaultTextFontSize',18);
set(sp2,'defaulttextinterpreter','latex');
set(sp2,'defaultAxesTickLabelInterpreter','latex');
set(sp2,'defaultLegendInterpreter','latex');
hold on; box on;
set(sp2,'fontsize',16);
set(gca,'Xticklabel',[])

% Plot yaw set-points
figure(6); colors = lines(4);
sp3 = subplot(3,1,3); grid on; box on; hold on;
% Lookup table bin to plot for example
wsi=1; % wind speed bin
tii=2; % turbulence intensity bin
for t=[1,3,4];
    % Number of samples in historical data
    n_full = squeeze(lookup_tables.n(wsi,tii,:,t));
    % Extract direction
    [dirw_full,indw_full] = sort(wrapTo180(lookup_tables.WndDir(n_full>=0)));
    % Extract yaw
    yaw_local = squeeze(lookup_tables.yaw(wsi,tii,n_full>=0,t));
    % Plot
    plot(dirw_full, yaw_local(indw_full),'o-', 'LineWidth', 1.5, 'color', colors(t,:), 'markerfacecolor', colors(t,:))
end
plot([-40,30],[0,0],'k-','LineWidth',1.75)
plot([-40,30],[0,0],'k-','LineWidth',1.75)
plot([xb(1),xb(1)],[-10000,10000],'k--','LineWidth',1.75);
plot([xb(2),xb(2)],[-10000,10000],'k--','LineWidth',1.75);
L=legend('$\mathrm{Turbine} \ 1$', ...
    '$\mathrm{Turbine} \ 2$', ...
    '$\mathrm{Turbine} \ 3$', ...
    'location','northwest','numcolumns',3);
set(L,'fontsize',14)
ylabel(strcat('$\gamma^* \ (^\circ)$'));
ylim([-30,30]); 
xlim([xb(1)-5, xb(2)+5]);
% Aesthetics
set(sp3,'DefaultAxesFontSize',18);
set(sp3,'DefaultTextFontSize',18);
set(sp3,'defaulttextinterpreter','latex');
set(sp3,'defaultAxesTickLabelInterpreter','latex');
set(sp3,'defaultLegendInterpreter','latex');
hold on; box on;
set(sp3,'fontsize',16);
pos3 = get(sp3, 'Position');
pos3_mod = pos3 - [-0.01 -0.05 0 0.05];
set(sp3, 'Position', pos3_mod )
xlabel('$\mathrm{Wind \ direction} \ \alpha \ (^\circ)$');