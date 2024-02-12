%% integration_segregation.m
%
%% Settings
function var_Barplot(bin_size,ttdays,regress_dist,outcome_option,option,conn,years)
close all;
paths;
ptList = {rns_config.patients.ID};
base_days = 90;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
if regress_dist
    suffix = '_regdist';
    figpath = fullfile(figpath,'regdist');
else
    suffix = '';
end
plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),suffix,'.mat'])).plasticity;
% Lead location x Connection type x Outcome loc x Freq
con_strings = {'intra','inter'};
depth_strings = {'Hippocampal','Neocortical'};
outcome_strings = ["Poor Responder","Good Responder"];
freq_strings = {'Theta','Alpha','Beta','Gamma'};
colors = ['g','b'];
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome_group}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
time = 90:bin_size:ttdays-1;
n_bin = length(time);
for i = 1:length(plasticity)
    plasticity(i).outcome_group = plasticity(i).outcome_group(outcome_option,:);
end
if ~exist(fullfile(figpath,'00_Network_trajectory',num2str(outcome_option)),'dir')
    mkdir(fullfile(figpath,'00_Network_trajectory',num2str(outcome_option)));
end
groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
ticks = groupCenters(2, 4, 1);
%% plot baseline si ratio
for y = 1:length(years)
    for conn = 1:2
        f = figure('Visible','Off');
        f.Position = [100,100,1200,600];
        for d = 1:2
            sub_data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
            [baseline_siratio_full,outcome] = getYearOutcome(sub_data,years{y},'outcome_group',option);
            baseline_siratio = cellfun(@(x,c) squeeze(mean(x(:,c==conn,:),2,'omitnan')),baseline_siratio_full,{sub_data.con_labels},'UniformOutput',false);
            baseline_siratio = horzcat(baseline_siratio{:})';
        
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
            if any(outcome == 1) && any(outcome == 2)
                for j = 1:4
                    p = ranksum(sr(:,j),pr(:,j));
                    text(ticks(j),ymax(2)-0.1,['p=',num2str(p,'%.3f')])
                end
            end
            ylabel(depth_strings{d})
            set(gca,'XTick',ticks,'XTickLabels',freq_strings)
            sgtitle(option)
        end
        saveas(f,fullfile(figpath,'00_Network_trajectory',num2str(outcome_option),[option,'_',num2str(years{y}),'_',con_strings{conn},'.png']))
    end
end
close all
end
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