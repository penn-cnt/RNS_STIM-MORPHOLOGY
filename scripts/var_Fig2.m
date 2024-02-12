% Plot Variance
%% Settings
function var_Fig2(bin_size,ttdays,regress_dist,outcome_option,option,years)
% variance of plv(%dplv) within time bins, plv_var/dplv_var
% variance across times bins, dplv_var_acrossbin/zplv_var_acrossbin
% variance in slope between time bins, plv_slope_var/dplv_slope_var
% variance in slope between events, dplv_slope_var_acrossbin/zplv_slope_var_acrossbin
% clear; close all;
paths;
ptList = {rns_config.patients.ID};
if regress_dist
    suffix = '_regdist';
    figpath = fullfile(figpath,'regdist');
else
    suffix = '';
end
plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),suffix,'.mat'])).plasticity;
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome_group}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
% Lead location x Connection type x Outcome loc x Freq
con_strings = ["Inside SOZ", "Between SOZ", "SOZ to Normal","Normal to Normal"];
con_strings2 = ["Within Lead", "Between Lead"];
depth_strings = {'Hippocampal','Neocortical'};
outcome_strings = ["Poor Responder","Good Responder"];
freq_strings = {'Theta','Alpha','Beta','Gamma'};
colors = ['g','b'];
% max_time_length = max(cellfun(@(x) size(x,2),dplv_traces),[],'all');
time = 90:bin_size:ttdays-1;
n_bin = length(time);
% uniform plv trace length

for i = 1:length(plasticity)
    plasticity(i).outcome_group = plasticity(i).outcome_group(outcome_option,:);
end
if ~exist(fullfile(figpath,'00_Network_trajectory',num2str(outcome_option),option),'dir')
    mkdir(fullfile(figpath,'00_Network_trajectory',num2str(outcome_option),option));
end
plasticity = paddata(plasticity,option,n_bin);
%% Figure 2 v2 within/between lead
ylims = [40,20];
addp = 5;
switch option
    case 'reorg_dplv'
        ylims = [40,20];
        addp = 5;
    case 'reorg_zplv'
        ylims = [2,1];
        addp = 0.1;
    case 'reorg_slope'
        ylims = [0.1,0.1];
        addp = 0.05;
end
for y = 1:length(years)
    for d = 1:2
        data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
        [dplv,outcome] = getYearOutcome(data,years{y},'outcome_group',option);
        f = figure('Visible','Off');
        f.Position = [100,100,1200,600];
        for i = 1:2 % conntype
            for j = 1:4 % freq
                subplot(2,4,sub2ind([4,2],j,i))
                hold on
                plv_traces = cellfun(@(x,conn) squeeze(x(:,conn == i,j)),dplv,{data.con_labels},'UniformOutput',false);
                for o = 1:2 % outcome
                    tmp = plv_traces(outcome == o);
                    traces = horzcat(tmp{:})';
                    if ~isempty(traces)
                        shadedErrorBar(time,mean(traces,1,'omitnan'), ...
                            std(traces,[],1,'omitnan')/sqrt(size(traces,1)), ...
                            "lineProps",['-',colors(o)],'transparent',1);
                        xlim([90,ttdays])
%                         ylim([-ylims(i),ylims(i)])
                    end
                end
                % FDA test and plot
                [fda_pv,sig_days] = FDA_test(plv_traces,outcome);
                text(ttdays-200,0,['p=',num2str(fda_pv,'%.3f')],'FontSize',12);
                if ~isempty(sig_days)
                    scatter(time(sig_days(:,1)),zeros(size(sig_days,1),1),[],'r','filled');
                end
                if i == 1
                    title(freq_strings{j},'FontWeight','bold')
                end
                if j == 1
                    ylabel(con_strings2(i),'FontWeight','bold')
                end
            end
        end
%         figure(d)
        legend(outcome_strings,'Position',[0.8,0.84,0.1,0.05]);
        sgtitle(depth_strings{d},'FontWeight','bold','FontSize',14)
        h = axes(f,'visible','off');
        h.XLabel.Visible='on';
        h.YLabel.Visible='on';
        ylabel(h,option,'Position',[-0.05,0.500000476837158,0], ...
                'FontWeight','bold','FontSize',14);
        xlabel(h,'Days','FontWeight','bold','FontSize',14);
        
        saveas(f,fullfile(figpath,'00_Network_trajectory',num2str(outcome_option),option,[depth_strings{d},'_',num2str(bin_size),'_',num2str(years{y}),'_lead.png']))
    end
    close all
end
end
