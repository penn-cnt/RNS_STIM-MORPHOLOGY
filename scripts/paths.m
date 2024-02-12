% Global settings of paths and necessary files
% global rootpath datapath rns_config
% Directories
rootpath = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(rootpath)); % add MORPHOLOGY-main folder
% default datapath
def_datapath = fullfile(rootpath,'user_data');
% load config file
try
    rns_config = jsondecode(fileread(fullfile(rootpath,"reference_data/config.JSON"))); 
catch
    error("No config file available.")
end
rns_toolbox_path = rns_config.paths.RNS_processing_toolbox;
assert(exist(rns_toolbox_path,'dir') == 7, 'Invalid RNS processing Toolbox path')
addpath(genpath(rns_toolbox_path))
if isfield(rns_config.paths,'RNS_PROCESSED_DATA_Folder')
    datapath = rns_config.paths.RNS_PROCESSED_DATA_Folder;
else
    datapath = def_datapath;
end
if ~exist(datapath, 'dir')
    mkdir(datapath);
end
figpath = fullfile(datapath,'figs');
if ~exist(figpath, 'dir')
    mkdir(figpath);
    mkdir(fullfile(figpath,'00_Network_trajectory'))
    mkdir(fullfile(figpath,'01_PLV_outcome'))
    mkdir(fullfile(figpath,'02_NTF_stim'))
    mkdir(fullfile(figpath,'03_Baseline_outcome'))
    mkdir(fullfile(figpath,'04_Baseline_plv'))
    mkdir(fullfile(figpath,'05_Baseline_stim'))
    mkdir(fullfile(figpath,'06_Int_seg'))
    mkdir(fullfile(figpath,'regdist'))
    mkdir(fullfile(figpath,'regdist','00_Network_trajectory'))
    mkdir(fullfile(figpath,'regdist','01_PLV_outcome'))
    mkdir(fullfile(figpath,'regdist','02_NTF_stim'))
    mkdir(fullfile(figpath,'regdist','03_Baseline_outcome'))
    mkdir(fullfile(figpath,'regdist','04_Baseline_plv'))
    mkdir(fullfile(figpath,'regdist','05_Baseline_stim'))
    mkdir(fullfile(figpath,'regdist','06_Int_seg'))
end
clear rootpath def_datapath rns_toolbox_path