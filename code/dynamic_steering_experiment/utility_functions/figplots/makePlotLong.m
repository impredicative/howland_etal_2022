function [ ] = makePlotBig(figNum)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
H = figure(figNum);
set(H,'color',[1 1 1],'InvertHardcopy','off','position',[50 50 1000 400]);
set(H,'DefaultAxesFontSize',18);
set(H,'DefaultTextFontSize',18);
set(H,'defaulttextinterpreter','latex');
set(H,'defaultAxesTickLabelInterpreter','latex');
set(H,'defaultLegendInterpreter','latex');
hold on; box on;
xlabel('$\left <u^{\prime 2} \right>$');
ylabel('$\left < u(x) u(x+r) \right >$');
set(gca,'fontsize',18);

end

