clear; close all; clc;
addpath("../reference_data/")
rns_config = jsondecode(fileread("reference_data/config.JSON")); %/Users/wojemann/Documents/CNT/RNS_STIM-MORPHOLOGY/
rootpath = rns_config.paths.RNS_project_code_path;
rns_toolbox_path = rns_config.paths.RNS_processing_toolbox;
addpath(genpath(rootpath))
addpath(genpath(rns_toolbox_path))
datapath = rns_config.paths.RNS_PROCESSED_DATA_Folder;
rns_config_path = fullfile(rns_toolbox_path, 'config.JSON');
patient_info = struct2table(load(fullfile(rns_config.paths.RNS_metadata_path,'patients_Penn.mat')).patients_Penn);

ptList = {rns_config.patients.ID};
nanzscore = @(x) (x - mean(x,2,"omitnan"))./std(x,[],[2,3],"omitnan");
fs = 250;
%% Trace Visualization (Figure 1a)

for pt = 1:length(ptList)

% Read Patient Data
ptID = ptList{pt};

analysis_windows = load([datapath,'/',ptID,'/stim_windows_',ptID,'.mat']).analysis_windows;
visit_selection_array = load([datapath,'/',ptID,'/visit_selection_array_',ptID,'.mat']).visit_selection_array;
ptime_trace = load([datapath,'/',ptID,'/posix_UTC_time_trace_',ptID,'.mat']).ptime_trace;

%% Separating Stimulations by stim
% Need to flip each channel so that the evoked response is negative
% [b,a] = butter(4,4/125,'high');
visits = unique(visit_selection_array);
all_visit_vals = {};
all_visit_times = {};
fprintf("%d unique events",length(visits))
for i_u = 1:length(visits)
    visit_stims = analysis_windows(visit_selection_array == visits(i_u));
    v_times = ptime_trace(visit_selection_array == visits(i_u));
    visit_vals = [];
    visit_times = [];
    for i_stim = 1:length(visit_stims)
        visit_data = visit_stims{i_stim}(750:1000,:);
%         visit_data = filtfilt(b,a,visit_data);
%         visit_vals(i_stim,:) = min(visit_stims{i_stim}(625:725,:),[],1);
        c = triu(corr(visit_data),1);
        try
            visit_vals(i_stim,:) = c(logical(c));
        catch
            visit_vals(i_stim,:) = zeros(1,6)*nan;
        end
        visit_times(i_stim) = v_times(i_stim);
    end
    all_visit_vals{i_u} = visit_vals;
    all_visit_times{i_u} = visit_times;
end

visit_lens = cellfun(@(x) size(x,1),all_visit_vals);
[~,i_sort] = sort(visit_lens,'descend');
sorted_visit_vals = all_visit_vals(i_sort);
sorted_visit_times = all_visit_times(i_sort);
sorted_lens = visit_lens(i_sort);
all_vals = zeros(length(visits),max(visit_lens),6); all_vals(:) = nan;
all_times = zeros(length(visits),max(visit_lens)); all_times(:) = nan;
for i = 1:length(visits)
    all_vals(i,1:sorted_lens(i),:) = sorted_visit_vals{i};
    all_times(i,1:sorted_lens(i)) = sorted_visit_times{i};
end
all_times = all_times-all_times(:,1);
% all_vals = nanzscore(all_vals);
%% Stim Response Evolution Plotting
figure(100 + pt)
corr_names = {"1 x 2","1 x 3","2 x 3","1 x 4", "2 x 4","3 x 4"};
for i_plot = 1:size(all_vals,3)
    subplot(3,2,i_plot)
    plot(all_times(1,:),all_vals(1,:,i_plot),'o')
    title(sprintf("Channels %s",corr_names{i_plot}))
end
sgtitle(ptID)
end
%% PCA
% all_stims = zeros(length(analysis_windows),length(analysis_windows{1}),4);
% for i_stim = 1:length(analysis_windows)
%     all_stims(i_stim,:,:) = analysis_windows{i_stim};
% end
% all_channels = [];
% for col = 1:4
%     all_channels = [all_channels; all_stims(:,:,col)];
% end
% figure
% coeff = pca(all_channels);
% scatter(coeff(1,:),coeff(2,:))
% % plot(all_stims(:,:,1)*coeff(:,3))
%% Trace Plotting
stim = 100;
[b,a] = butter(4,4/125,'high');
figure(pt)
ax = gca;
data = analysis_windows{stim};
filt_data = filtfilt(b,a,data);
plot(filt_data+(4:-1:1)*500)
xlabel('Time (Seconds)','FontSize',18)

ylabel('Signal Recording','Fontsize',18)
yticks((1:4)*500);
yticklabels({'Channel 4','Channel 3','Channel 2','Channel 1'})
f = gcf;
ax.FontSize = 16;
title(['Interictal Stimulation: ' ptID],'FontSize',24)
% exportgraphics(gcf,'stim_trace.pdf')