%% baseline_stim.m
baseline_stim_pred = load(fullfile(datapath,'NTF_stim_pred.mat')).NTF_stim_pred;
for pt = 1:length(baseline_stim_pred)
    ptID = baseline_stim_pred(pt).ID;
    p = find(strcmp(ptID,{plasticity.ptID}));
    baseline_stim_pred(pt).baseline_plvs = plasticity(p).baseline_plvs(1,:);
end
%% Fig.4C
freq_strings = {'Theta','Alpha','Beta','Gamma'};
windows = [0, 1, 2, 3, 4, 7, 14, 30, 90, 180, 360, 720];
outcome = cellfun(@(x) x(1), {baseline_stim_pred.outcome});
R2_combined = vertcat(baseline_stim_pred.R2);
bplv_combined = vertcat(baseline_stim_pred.baseline_plvs);
% R2_all = [];
% params_all = {};
% preds_all = {};
% trues_all = {};
% R2_all_null = {};
% pR2_all = [];
for w = 1:length(windows)
    mdl = fitglm(bplv_combined,R2_combined(:,w));
    p = mdl.Coefficients{:,"pValue"};
    R2 = mdl.Rsquared.Ordinary;
    mdl = fitglm(bplv_combined(outcome == 1,:),R2_combined(outcome == 1,w));
    p_PR = mdl.Coefficients{:,"pValue"};
    R2_PR = mdl.Rsquared.Ordinary;
    mdl = fitglm(bplv_combined(outcome == 2,:),R2_combined(outcome == 2,w));
    p_SR = mdl.Coefficients{:,"pValue"};
    R2_SR = mdl.Rsquared.Ordinary;
    f = figure('Position',[100,100,1200,300]);
    for freq = 1:4
        subplot(1,4,freq)
        scatter(bplv_combined(outcome == 1,freq),R2_combined(outcome == 1,w),[],'g','filled');
        hold on
        scatter(bplv_combined(outcome == 2,freq),R2_combined(outcome == 2,w),[],'b','filled');
        xlim([-0.05,1.05])
        ylim([-0.05,1.05])
        title(freq_strings(freq),'FontWeight','bold','FontSize',12)
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
    saveas(f,fullfile(datapath,'figs','baseline_stim',['window_',num2str(windows(w)),'.png']))
end
