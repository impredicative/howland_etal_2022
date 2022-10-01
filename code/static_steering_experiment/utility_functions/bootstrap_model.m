function [stats] = bootstrap_model(data)
% Bootstrap model output
data.inds = 1:length(data.p);
stats={};

% Bootstrapping
cut = (1-data.conf)/2;
ind_low = ceil(data.Boots*cut); ind_high = round(data.Boots*(1-cut));
Prm_all_local = zeros(data.Boots,1);
for i=1:data.Boots;
    inds_local = randsample(data.inds,length(data.inds),true);
    Prm_all_local(i) = nanmean(data.p(inds_local));
end
Prm_all_local = sort(Prm_all_local,1);
stats.P_CI(1) = mean(Prm_all_local) - Prm_all_local(ind_low);
stats.P_CI(2) = Prm_all_local(ind_high) - mean(Prm_all_local);
stats.P_mean = mean(Prm_all_local);

end