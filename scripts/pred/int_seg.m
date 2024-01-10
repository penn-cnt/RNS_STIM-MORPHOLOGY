%% integration_segregation.m
% 
%% Settings
clear; close all;
paths;
ptList = {rns_config.patients.ID};
bin_size = 90; % subject to change
ttdays = 1050;
plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),'.mat'])).plasticity;
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome_group}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
% Lead location x Connection type x Outcome loc x Freq
con_strings = ["Inside SOZ", "Between SOZ", "SOZ to Normal","Normal to Normal"];
depth_strings = {'Hippocampal','Neocortical'};
outcome_strings = ["Poor Responder","Good Responder"];
freq_strings = {'Theta','Alpha','Beta','Gamma'};
colors = ['g','b'];
% max_time_length = max(cellfun(@(x) size(x,2),dplv_traces),[],'all');
time = 90:bin_size:ttdays-1;
n_bin = length(time);
for n = 1:length(plasticity)
    try
        plasticity(n).reorg_dplv = cellfun(@(x) x(1:n_bin,:),plasticity(n).reorg_dplv,'UniformOutput',false);
    catch
        plasticity(n).reorg_dplv = cellfun(@(x) padarray(x,n_bin-size(x,1),'post'),plasticity(n).reorg_dplv,'UniformOutput',false);
    end
    plasticity(n).seg_int = [];
    for f = 1:4
        tmp = mean(diff(plasticity(n).reorg_dplv{1,f},1,1),'all','omitnan');
        tmp2 =  mean(diff(plasticity(n).reorg_dplv{2,f},1,1),'all','omitnan');
        plasticity(n).seg_int = [plasticity(n).seg_int, tmp2 - tmp];
    end
end
groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
ticks = groupCenters(2, 4, 1);
%% Figure 2 v2
years = {1,2,3,'end'};
for y = 1:length(years)
    if isnumeric(years{y})
        sub_data = plasticity(cellfun(@(x) length(x) >= years{y}, {plasticity.outcome_group}));
        outcome = cellfun(@(x) x(years{y}), {sub_data.outcome_group});
        year = num2str(years{y});
    elseif strcmp(years{y},'end')
        sub_data = plasticity;
        outcome = cellfun(@(x) x(end), {sub_data.outcome_group});
        year = 'end';
    end
    f = figure(y);
    f.Position = [100,100,1200,600];
    for d = 1:2
        subplot(2,1,d)
        hold on
        ind = cellfun(@(x) x==d,{sub_data.depth});
        data = vertcat(sub_data(ind).seg_int); 
        o = outcome(ind);
        box_data = {};
        sr = data(o == 1,:);
        pr = data(o == 2,:);
        sr = padarray(sr,[max(size(sr,1),size(pr,1))-size(sr,1),0],nan,'post');
        pr = padarray(pr,[max(size(sr,1),size(pr,1))-size(pr,1),0],nan,'post');
        box_data{1} = pr;
        box_data{2} = sr;
        b = boxplotGroup(box_data,'primaryLabels',repmat({''}, 2, 1), ...
            'Colors',[0,1,0;0,0,1],'GroupType','betweenGroups', ...
            'PlotStyle','traditional','BoxStyle','outline', ...
            'Symbol','o','Widths',0.7);
        ymax = get(gca, 'YLim');
        for j = 1:4
            p = ranksum(sr(:,j),pr(:,j));
            text(ticks(j),ymax(2)-2,['p=',num2str(p,'%.3f')])
        end
        ylabel(depth_strings{d})
        set(gca,'XTick',ticks,'XTickLabels',freq_strings)
        sgtitle('Int-seg Difference')
    end
    saveas(f,fullfile(figpath,'06_Int_seg',[year,'.png']))
end
close all