%% dplv_traces.m
% This script replicates figure 2 of Ankit's paper.
%% Settings
clear; close all;
paths;
ptList = {rns_config.patients.ID};
bin_size = 90; % subject to change
ttdays = 1050; % total days to plot traces
plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),'.mat'])).plasticity;
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome_group}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
% Lead location x Connection type x Outcome loc x Freq
con_strings = ["Inside SOZ", "Between SOZ", "SOZ to Normal","Normal to Normal"];
con_strings2 = ["Within Lead", "Between Lead"];
depth_strings = {'Hippocampal','Neocortical'};
outcome_strings = ["Poor Responder","Good Responder"];
freq_strings = {'Theta','Alpha','Gamma','Beta'};
colors = ['g','b'];
% max_time_length = max(cellfun(@(x) size(x,2),dplv_traces),[],'all');
time = 90:bin_size:ttdays-1;
n_bin = length(time);
% uniform plv trace length
for n = 1:length(plasticity)
    try
        plasticity(n).reorg_dplv = cellfun(@(x) x(1:n_bin,:),plasticity(n).reorg_dplv,'UniformOutput',false);
    catch
        plasticity(n).reorg_dplv = cellfun(@(x) padarray(x,n_bin-size(x,1),nan,'post'),plasticity(n).reorg_dplv,'UniformOutput',false);
    end
end
years = {1,2,3,'end'};
%% Figure 2 v1 soz-soz 
ylims = [40,20];
for y = 1:length(years)
    if isnumeric(years{y})
        all_year = num2cell(years{y} * ones(1,length(plasticity)));
        year = num2str(years{y});
    elseif strcmp(years{y},'end')
        year = 'end';
        all_year = num2cell(cellfun(@(x) length(x),{plasticity.outcome}));
    end
    for d = 1:2
        f = figure(d);
        f.Position = [100,100,1200,600];
        for i = 1:2 % conntype
            for j = 1:4 % freq
                subplot(2,4,sub2ind([4,2],j,i))
                hold on
                ind = cellfun(@(x) x==d,{plasticity.depth}) & ...
                      cellfun(@(x) length(x) >= all_year{1},{plasticity.outcome_group});
                data = plasticity(ind);
                tmp_year = all_year(ind);
                for n = 1:length(data)
                    data(n).outcome = data(n).outcome_group(tmp_year{n});
                    data(n).reorg_dplv = data(n).reorg_dplv{i+2,j};
                end
                for o = 1:2 % outcome
                    tmp = data(cellfun(@(x,y) x(y)==o,{data.outcome_group},tmp_year));
                    traces = horzcat(tmp.reorg_dplv)';
                    if ~isempty(traces)
                        traces = traces(:,1:length(time));
                        traces = traces(exclude_outlier(traces,1.5),:);
                        shadedErrorBar(time,mean(traces,1,'omitnan'), ...
                            std(traces,[],1,'omitnan')/sqrt(size(traces,1)), ...
                            "lineProps",['-',colors(o)],'transparent',1);
                        ylabel(con_strings2(i))
                        xlim([90,ttdays])
                        ylim([-ylims(i),ylims(i)])
                    end
                end
                % FDA test and plot
                outcome = [data.outcome];
                data = {data.reorg_dplv};
                [fda_pv,sig_days] = FDA_test(data,outcome);
                text(ttdays-200,-ylims(i)+5,['p=',num2str(fda_pv,'%.3f')],'FontSize',12);
                if ~isempty(sig_days)
                    scatter(time(sig_days(:,1)),zeros(size(sig_days,1),1),[],'r','filled');
                end
            end
        end
        figure(d)
        sgtitle(depth_strings{d})
        legend(outcome_strings,'Position',[0.8,0.87,0.1,0.05]);
        saveas(f,fullfile(figpath,'00_Network_trajectory',[depth_strings{d},'_',num2str(bin_size),'_',year,'_soz.png']))
    end
    close all
end
%% Figure 2 v2 within/between lead
ylims = [40,20];
for y = 1:length(years)
    for d = 1:2
        data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
        [dplv,outcome] = getYearOutcome(data,years{y},'outcome_group','reorg_dplv');
        f = figure(d);
        f.Position = [100,100,1200,600];
        for i = 1:2 % conntype
            for j = 1:4 % freq
                subplot(2,4,sub2ind([4,2],j,i))
                hold on
                plv_traces = cellfun(@(x) x{i,j},dplv,'UniformOutput',false);
                for o = 1:2 % outcome
                    tmp = plv_traces(outcome == o);
                    traces = horzcat(tmp{:})';
                    if ~isempty(traces)
                        shadedErrorBar(time,mean(traces,1,'omitnan'), ...
                            std(traces,[],1,'omitnan')/sqrt(size(traces,1)), ...
                            "lineProps",['-',colors(o)],'transparent',1);
                        xlim([90,ttdays])
                        ylim([-ylims(i),ylims(i)])
                    end
                end
                % FDA test and plot
                [fda_pv,sig_days] = FDA_test(plv_traces,outcome);
                text(ttdays-200,-ylims(i)+5,['p=',num2str(fda_pv,'%.3f')],'FontSize',12);
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
        figure(d)
        legend(outcome_strings,'Position',[0.8,0.84,0.1,0.05]);
        sgtitle(depth_strings{d},'FontWeight','bold','FontSize',14)
        h = axes(f,'visible','off');
        h.XLabel.Visible='on';
        h.YLabel.Visible='on';
        ylabel(h,'% dPLV','Position',[-0.05,0.500000476837158,0], ...
                'FontWeight','bold','FontSize',14);
        xlabel(h,'Days','FontWeight','bold','FontSize',14);
        
        saveas(f,fullfile(figpath,'00_Network_trajectory',[depth_strings{d},'_',num2str(bin_size),'_',num2str(years{y}),'_lead.png']))
    end
    close all
end
%%
figure(3);
for i = 1:4 % conntype
    for j = 1:4 % freq
        subplot(4,4,sub2ind([4,4],j,i))
        hold on
        for o = 1:2 % outcome
            trace = [dplv_traces{1,i,o,j};dplv_traces{2,i,o,j}];
            if ~isempty(trace)
                shadedErrorBar(time,median(trace,1,'omitnan'), ...
                    std(trace,[],1,'omitnan')/sqrt(size(trace,1)), ...
                    "lineProps",['-',colors(o)],'transparent',1);
                ylabel(con_strings(i))
                xlim([90,800])
                ylim([-0.5,0.5])
            end
        end
    end
end

% if plot_figure
%     % get plv instead of dplv time trace
%
%
% end
% if plot_figure
%     figure(pt)
%     hold on
%     for i = 1:size(band_limited_plvs,2) % each connection
%         for j = 1:size(dplv,3) % each frequency
%             subplot(7,1,i)
%             hold on
%             plot(binned_times(logical(binned_times)),binned_dplvs(:,i,j))
%             xlim([90,max(binned_times)])
%             mdl = fitglm(test(logical(test)),binned_stims(logical(binned_stims)),"y~x1","Distribution","poisson");
%             title(num2str(mdl.Rsquared.Ordinary))
%         end
%     end
%     subplot(7,1,7)
%     plot(binned_times(logical(binned_times)),binned_counts(logical(binned_counts)))
%     xlim([90,max(binned_times)])
%     sgtitle(ptID)
%     %% Network NMF Trajectories
%     [W,H,Q] = betaNTF(dplv,3);
%     % figure
%     % plot(time_trace(~baseline_mask),W(:,3))
%     mdl = fitglm(W(:,3),stim_trace(~baseline_mask)','y ~ x1','Distribution',"poisson");
%     R2 = mdl.Rsquared.Ordinary;
% end
%

