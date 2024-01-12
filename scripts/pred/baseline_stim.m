%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;
plasticity = load(fullfile(datapath,"plasticity_90.mat")).plasticity;
%% baseline_stim.m
baseline_stim_pred = load(fullfile(datapath,'NTF_stim_pred.mat')).NTF_stim_pred;
for pt = 1:length(baseline_stim_pred)
    ptID = baseline_stim_pred(pt).ID;
    p = find(strcmp(ptID,{plasticity.ptID}));
    baseline_stim_pred(pt).baseline_plvs = plasticity(p).baseline_plvs;
end
%% Fig.4C
freq_strings = {'Theta','Alpha','Beta','Gamma'};
freq_strings2 = {'Fac1','Fac2','Fac3'};
windows = [0, 1, 2, 3, 4, 7, 14, 30, 90, 180, 360, 720];
years = {1,2,3,'end'};
data = baseline_stim_pred(~cellfun('isempty',{baseline_stim_pred.outcome}));
for y = 1:length(years)
    [all_data,outcome] = getYearOutcome(data,years{y},'outcome',{'R2','baseline_plvs','base_fac_time'});
    R2_combined = vertcat(all_data{1}{:});
    bplv_intra_combined = cellfun(@(x) mean(x(1:2,:),'omitnan'), all_data{2},'UniformOutput',false);
    bplv_intra_combined = vertcat(bplv_intra_combined{:});
    bplv_inter_combined = cellfun(@(x) mean(x(3:6,:),'omitnan'), all_data{2},'UniformOutput',false);
    bplv_inter_combined = vertcat(bplv_inter_combined{:});
    bfc_combined = vertcat(all_data{3}{:});
    all_data = {bplv_intra_combined,bplv_inter_combined,bfc_combined};
    label = {'intra','inter','bfc'};
    for n = 1:3
        tmp = all_data{n};
        for w = 1:length(windows)
            mdl = fitglm(tmp,R2_combined(:,w));
            p = mdl.Coefficients{:,"pValue"};
            R2 = mdl.Rsquared.Ordinary;
            mdl = fitglm(tmp(outcome == 1,:),R2_combined(outcome == 1,w));
            p_PR = mdl.Coefficients{:,"pValue"};
            R2_PR = mdl.Rsquared.Ordinary;
            mdl = fitglm(tmp(outcome == 2,:),R2_combined(outcome == 2,w));
            p_SR = mdl.Coefficients{:,"pValue"};
            R2_SR = mdl.Rsquared.Ordinary;
            f = figure('Position',[100,100,1200,300]);
            for freq = 1:size(tmp,2)
                subplot(1,size(tmp,2),freq)
                scatter(tmp(outcome == 1,freq),R2_combined(outcome == 1,w),[],'g','filled');
                hold on
                scatter(tmp(outcome == 2,freq),R2_combined(outcome == 2,w),[],'b','filled');
                xlim([-0.05,1.05])
                ylim([-0.05,1.05])
                if size(tmp,2) == 4
                    title(freq_strings(freq),'FontWeight','bold','FontSize',12)
                else 
                    title(freq_strings2(freq),'FontWeight','bold','FontSize',12)
                end
                xlabel('Baseline PLV','FontWeight','bold','FontSize',12)
                if freq == 1
                    ylabel('NTF-Stim R2','FontWeight','bold','FontSize',12)
                end
                text(0.8,0.8,{['p=',num2str(p(freq+1),'%.3f')], ...
                    ['p_{SR}=',num2str(p_SR(freq+1),'%.3f')], ...
                    ['p_{PR}=',num2str(p_PR(freq+1),'%.3f')]},'FontSize',12)
            end
            sgtitle(['Window Size ',num2str(windows(w)), ...
                ' R^2 = ',num2str(R2,'%.2f'),' | ' , ...
                num2str(R2_SR,'%.2f'),' | ' ,...
                num2str(R2_PR,'%.2f')],'FontWeight','bold','FontSize',14)
            saveas(f,fullfile(figpath,'05_Baseline_stim',label{n},['year_',num2str(years{y}),'_window_',num2str(windows(w)),'.png']))
        end
        close all
    end
end

