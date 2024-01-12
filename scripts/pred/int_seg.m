%% integration_segregation.m
%
%% Settings
clear; close all;
paths;
ptList = {rns_config.patients.ID};
base_days = 90;
bin_size = 90; % subject to change
ttdays = 1050;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),'.mat'])).plasticity;
% Lead location x Connection type x Outcome loc x Freq
con_strings = ["Inside SOZ", "Between SOZ", "SOZ to Normal","Normal to Normal"];
depth_strings = {'Hippocampal','Neocortical'};
outcome_strings = ["Poor Responder","Good Responder"];
freq_strings = {'Theta','Alpha','Beta','Gamma'};
colors = ['g','b'];
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome_group}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
time = 90:bin_size:ttdays-1;
n_bin = length(time);
years = {1,2,3,'end'};
for n = 1:length(plasticity)
    try
        plasticity(n).reorg_dplv = cellfun(@(x) x(1:n_bin,:),plasticity(n).reorg_dplv,'UniformOutput',false);
        plasticity(n).si_ratio = plasticity(n).si_ratio(1:n_bin,:);
        plasticity(n).int_seg = plasticity(n).int_seg(1:n_bin,:);
    catch
        plasticity(n).reorg_dplv = cellfun(@(x) padarray(x,n_bin-size(x,1),nan,'post'),plasticity(n).reorg_dplv,'UniformOutput',false);
        plasticity(n).si_ratio = padarray(plasticity(n).si_ratio,n_bin-size(plasticity(n).si_ratio,1),nan,'post');
        plasticity(n).int_seg = padarray(plasticity(n).int_seg,n_bin-size(plasticity(n).int_seg,1),nan,'post');
    end
end
groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
ticks = groupCenters(2, 4, 1);
%% plot si ratio
ylims = [40,20];
for y = 1:length(years)
    f = figure('Position',[100,100,1200,600]);
    for d = 1:2
        data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
        [si_ratio,outcome] = getYearOutcome(data,years{y},'outcome_group','si_ratio');
        for j = 1:4 % freq
            subplot(2,4,sub2ind([4,2],j,d))
            hold on
            data_freq = cellfun(@(x) x(:,j),si_ratio,'UniformOutput',false);
            for o = 1:2 % outcome
                tmp = data_freq(outcome == o);
                traces = horzcat(tmp{:})';
                if ~isempty(traces)
                    shadedErrorBar(time,mean(traces,1,'omitnan'), ...
                        std(traces,[],1,'omitnan')/sqrt(size(traces,1)), ...
                        "lineProps",['-',colors(o)],'transparent',1);
                    xlim([90,ttdays])
                    ylim([0,1])
                end 
            end
            % FDA test and plot
            [fda_pv,sig_days] = FDA_test(data_freq,outcome);
            text(ttdays-200,0.2,['p=',num2str(fda_pv,'%.3f')],'FontSize',12);
            if ~isempty(sig_days)
                scatter(time(sig_days(:,1)),0.1*ones(size(sig_days,1),1),[],'r','filled');
            end
            if d == 1
                title(freq_strings{j},'FontWeight','bold')
            end
            if j == 1
                ylabel(depth_strings(d),'FontWeight','bold')
            end
        end
    end
    legend(outcome_strings,'Position',[0.8,0.84,0.1,0.05]);
    h = axes(f,'visible','off');
    h.XLabel.Visible='on';
    h.YLabel.Visible='on';
    ylabel(h,'Int-Seg Ratio','Position',[-0.05,0.500000476837158,0], ...
        'FontWeight','bold','FontSize',14);
    xlabel(h,'Days','FontWeight','bold','FontSize',14);
    saveas(f,fullfile(figpath,'06_Int_seg',['si_ratio_',num2str(years{y}),'.png']))
end
close all
%% plot int-seg dPLV traces
ylims = [40,20];
for y = 1:length(years)
    f = figure('Position',[100,100,1200,600]);
    for d = 1:2
        data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
        [seg_int,outcome] = getYearOutcome(data,years{y},'outcome_group','int_seg');
        for j = 1:4 % freq
            subplot(2,4,sub2ind([4,2],j,d))
            hold on
            data_freq = cellfun(@(x) x(:,j),seg_int,'UniformOutput',false);
            for o = 1:2 % outcome
                tmp = data_freq(outcome == o);
                traces = horzcat(tmp{:})';
                if ~isempty(traces)
                    shadedErrorBar(time,mean(traces,1,'omitnan'), ...
                        std(traces,[],1,'omitnan')/sqrt(size(traces,1)), ...
                        "lineProps",['-',colors(o)],'transparent',1);
                    xlim([90,ttdays])
                    ylim([-60,60])
                end
            end
            % FDA test and plot
            [fda_pv,sig_days] = FDA_test(data_freq,outcome);
            text(ttdays-200,-55,['p=',num2str(fda_pv,'%.3f')],'FontSize',12);
            if ~isempty(sig_days)
                scatter(time(sig_days(:,1)),zeros(size(sig_days,1),1),[],'r','filled');
            end
            if d == 1
                title(freq_strings{j},'FontWeight','bold')
            end
            if j == 1
                ylabel(depth_strings(d),'FontWeight','bold')
            end
        end
    end
    legend(outcome_strings,'Position',[0.8,0.84,0.1,0.05]);
    h = axes(f,'visible','off');
    h.XLabel.Visible='on';
    h.YLabel.Visible='on';
    ylabel(h,'Int-Seg Difference','Position',[-0.05,0.500000476837158,0], ...
        'FontWeight','bold','FontSize',14);
    xlabel(h,'Days','FontWeight','bold','FontSize',14);
    saveas(f,fullfile(figpath,'06_Int_seg',['traces_',num2str(years{y}),'.png']))
end
close all
%% plot baseline si ratio
for y = 1:length(years)
    f = figure(y);
    f.Position = [100,100,1200,600];
    for d = 1:2
        sub_data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
        [baseline_siratio,outcome] = getYearOutcome(sub_data,years{y},'outcome_group','baseline_siratio');
        baseline_siratio = vertcat(baseline_siratio{:});
        subplot(2,1,d)
        hold on
        box_data = {};
        sr = baseline_siratio(outcome == 2,:);
        pr = baseline_siratio(outcome == 1,:);
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
            text(ticks(j),ymax(2)-0.1,['p=',num2str(p,'%.3f')])
        end
        ylabel(depth_strings{d})
        set(gca,'XTick',ticks,'XTickLabels',freq_strings)
        sgtitle('Int-seg Difference')
    end
    saveas(f,fullfile(figpath,'06_Int_seg',['baseline_siratio_',num2str(years{y}),'.png']))
end
close all
%% plot diff dPLV
% for y = 1:length(years)
%     f = figure(y);
%     f.Position = [100,100,1200,600];
%     for d = 1:2
%         sub_data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
%         [baseline_siratio,outcome] = getYearOutcome(sub_data,years{y},'outcome_group','seg_int');
%         baseline_siratio = vertcat(baseline_siratio{:});
%         subplot(2,1,d)
%         hold on
%         box_data = {};
%         sr = baseline_siratio(outcome == 2,:);
%         pr = baseline_siratio(outcome == 1,:);
%         sr = padarray(sr,[max(size(sr,1),size(pr,1))-size(sr,1),0],nan,'post');
%         pr = padarray(pr,[max(size(sr,1),size(pr,1))-size(pr,1),0],nan,'post');
%         box_data{1} = pr;
%         box_data{2} = sr;
%         b = boxplotGroup(box_data,'primaryLabels',repmat({''}, 2, 1), ...
%             'Colors',[0,1,0;0,0,1],'GroupType','betweenGroups', ...
%             'PlotStyle','traditional','BoxStyle','outline', ...
%             'Symbol','o','Widths',0.7);
%         ymax = get(gca, 'YLim');
%         for j = 1:4
%             p = ranksum(sr(:,j),pr(:,j));
%             text(ticks(j),ymax(2)-2,['p=',num2str(p,'%.3f')])
%         end
%         ylabel(depth_strings{d})
%         set(gca,'XTick',ticks,'XTickLabels',freq_strings)
%         sgtitle('Int-seg Difference')
%     end
%     saveas(f,fullfile(figpath,'06_Int_seg',[num2str(years{y}),'.png']))
% end
% close all