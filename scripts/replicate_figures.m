%% Plotting for replicating Ankit's Fig 2
if plot_figure
    % get plv instead of dplv time trace


end
if plot_figure
    figure(pt)
    hold on
    for i = 1:size(band_limited_plvs,2) % each connection
        for j = 1:size(dplv,3) % each frequency
            subplot(7,1,i)
            hold on
            plot(binned_times(logical(binned_times)),binned_dplvs(:,i,j))
            xlim([90,max(binned_times)])
            mdl = fitglm(test(logical(test)),binned_stims(logical(binned_stims)),"y~x1","Distribution","poisson");
            title(num2str(mdl.Rsquared.Ordinary))
        end
    end
    subplot(7,1,7)
    plot(binned_times(logical(binned_times)),binned_counts(logical(binned_counts)))
    xlim([90,max(binned_times)])
    sgtitle(ptID)
    %% Network NMF Trajectories
    [W,H,Q] = betaNTF(dplv,3);
    % figure
    % plot(time_trace(~baseline_mask),W(:,3))
    mdl = fitglm(W(:,3),stim_trace(~baseline_mask)','y ~ x1','Distribution',"poisson");
    R2 = mdl.Rsquared.Ordinary;
end




%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;
bin_size = 45; % subject to change

%% Plotting
plv_traces = load(fullfile(datapath,"plv_traces.mat")).plv_traces;
dplv_traces = load(fullfile(datapath,"plv_traces.mat")).dplv_traces;
% Lead location x Connection type x Outcome loc x Freq
con_strings = ["inside soz", "between soz", "soz to normal","normal to normal"];
depth_strings = ["Hippocampal","Neocortical"];
colors = ['r','b'];
max_time_length = size(plv_traces{1,2,1,1},2);
time = bin_size:bin_size:bin_size*max_time_length;
time = time + 90;
%%
for d = 1:2
    figure(d);
    for i = 1:4 % conntype
        for j = 1:4 % freq
            subplot(4,4,sub2ind([4,4],j,i))
            hold on
            for o = 1:2 % outcome
                if ~isempty(plv_traces{d,i,o,j})
                    shadedErrorBar(time,median(dplv_traces{d,i,o,j},1,'omitnan'), ...
                        std(plv_traces{d,i,o,j},[],1,'omitnan'), ...
                        "lineProps",['-',colors(o)],'transparent',1);
                    ylabel(con_strings(i))
                    xlim([90,900])
                    ylim([-0.5,0.5])
                end
            end
        end
    end
    figure(d)
    sgtitle(depth_strings(d))
end
