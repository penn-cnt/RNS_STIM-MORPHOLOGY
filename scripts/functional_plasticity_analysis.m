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
%% Connectivity Trajectories

for pt = 1:length(ptList)
%% Read Patient Data
ptID = ptList{pt};
disp(ptID)
time_trace = load([datapath,'/',ptID,'/UTC_time_trace_',ptID,'.mat']).time_trace;
ptime_trace = load([datapath,'/',ptID,'/posix_UTC_time_trace_',ptID,'.mat']).ptime_trace;
all_plvs = load([datapath,'/',ptID,'/plvs_',ptID,'.mat']).all_plvs;


[~,~,~, histT] = loadRNSptData(ptID, rns_config);
stim_data = histT(:,["UTCStartTime","EpisodeStarts"]);
cum_stims = movingCumulativeSum(stim_data{:,"EpisodeStarts"},90*24);
[~,stim_idxs] = pdist2(posixtime(stim_data{:,"UTCStartTime"}),ptime_trace,"euclidean",'Smallest',1);
stim_trace = cum_stims(stim_idxs);

freq_array = 4:100;
freqs = [freq_array(1:end-1)',freq_array(2:end)'];

%% Band Limited Network Trajectories
% Need to examine the lowest time resolution that we can find within this
% patient. Eventually may need to find lowest across patients.
% plot(diff(time_trace))
baseline_period = time_trace(1) + days(90);
baseline_mask = time_trace < baseline_period;
implant_time = time_trace-time_trace(1);
implant_time = days(implant_time(~baseline_mask));

freq_bands_mask = {(freqs(:,1) > 0 & freqs(:,2) <= 8),...
    (freqs(:,1) > 8) & (freqs(:,2) <= 15),...
    (freqs(:,1) > 15) & (freqs(:,2) <= 30),...
    (freqs(:,1) > 30) & (freqs(:,2) <= 100)};
band_limited_plvs = zeros(size(all_plvs,1),size(all_plvs,2),length(freq_bands_mask));
for mask = 1:length(freq_bands_mask)
    band_limited_plvs(:,:,mask) = mean(all_plvs(:,:,mask),3);
end
baseline_values = mean(band_limited_plvs(baseline_mask,:,:),1);

dplv = (band_limited_plvs(~baseline_mask,:,:)-baseline_values)./baseline_values;

bin_size = 60;
bin_indices = floor(implant_time/bin_size);
binned_counts = accumarray(bin_indices,1,[],@sum);
binned_times = accumarray(bin_indices,implant_time,[],@max);
binned_stims = accumarray(bin_indices,stim_trace(~baseline_mask));
binned_dplvs = zeros(length(unique(bin_indices)),size(dplv,2),size(dplv,3));

figure(pt)
hold on
for i = 1:size(band_limited_plvs,2)
    for j = 1:size(dplv,3)
        try
        test = accumarray(bin_indices,dplv(:,i,j),[],@mean);
        binned_dplvs(:,i,j) = test(logical(test));
        subplot(6,1,i)
        hold on
        plot(binned_times(logical(binned_times)),test(logical(test)))
        xlim([90,810])
        mdl = fitglm(test(logical(test)),binned_stims(logical(test)),"y~x1","Distribution","poisson");
        title(num2str(mdl.Rsquared.Ordinary))
        catch
            disp("likely nan's in data")
        end
    end
end
continue
%% Network NMF Trajectories
[W,H,Q] = betaNTF(dplv,3);
%%
% figure
% plot(time_trace(~baseline_mask),W(:,3))
mdl = fitglm(W(:,3),stim_trace(~baseline_mask)','y ~ x1','Distribution',"poisson")
R2 = mdl.Rsquared.Ordinary

end
%% Functions
function movingSum = movingCumulativeSum(inputArray, windowSize)
    % Check if the window size is valid
    if windowSize <= 0 || windowSize > length(inputArray)
        error('Window size must be a positive integer smaller than or equal to the length of the input array.');
    end
    
    % Initialize the output array to store the moving cumulative sum
    movingSum = zeros(1, length(inputArray));
    
    % Calculate the cumulative sum for the first window
    movingSum(1:windowSize) = cumsum(inputArray(1:windowSize),"omitnan");
    
    % Iterate through the rest of the array to update the moving sum
    for i = windowSize+1:length(inputArray)
        % Update the moving sum by adding the current element and subtracting the element that is now out of the window
        movingSum(i) = sum(inputArray((i-windowSize+1):i),"omitnan");
%         movingSum(i) = movingSum(i-1) + inputArray(i) - inputArray(i - windowSize);
        if isnan(movingSum(i))
            disp('x')
        end
    end
end

function scaledData = scaleDataToMinus1To1(data)
    % Find the minimum and maximum values in the data
    minValue = min(data);
    maxValue = max(data);
    
    % Scale the data to the range [-1, 1]
    scaledData = ((data - minValue) * 2) ./ (maxValue - minValue) - 1;
end