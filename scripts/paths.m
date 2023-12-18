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