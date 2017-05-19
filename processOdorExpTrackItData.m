% A wrapper script file to invoke all the track it related functions in the
% right order.

%% Defaults
force_rewrite = 0;

inputDir_pref = '/media/dinesh/Odor Conc Exp Data/';
outputDir_pref = '/media/dinesh/Odor Conc Exp Data/';
moveDir_pref = '/media/dinesh/Odor Conc Exp Data/';

%% Get the experiment folder, output folder and backup folder
expDir = uigetdir(inputDir_pref,'Experiment Folder (Trackit save folder)');
outDir = uigetdir(outputDir_pref,'Output Folder (Main Data folder)');
backupDir = uigetdir(moveDir_pref,'Backup Folder (To move data)');

%% Process trackit data
raw_data_matfile = copyTrackitExperimentData(expDir,outDir,backupDir);
all_traj_matfile = extractTrajectories(raw_data_matfile, force_rewrite);
sorted_matfile = extractFirstLandedTrajectories(all_traj_matfile, force_rewrite);

close all;
