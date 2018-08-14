% A wrapper script file to invoke all the track it related functions in the
% right order.
% 
% Dinesh Natesan 


%% Defaults
force_rewrite = 1;
backup_on = 0;

%% Experiment folders (hardcoded)
expDir = '/media/dinesh/M2SSD/Odor Concentration Data/Airflow';
outDir = '/media/dinesh/M2SSD/Odor Concentration Analysis/Airflow/Output';
backupDir = '/media/dinesh/M2SSD/Odor Concentration Analysis/Airflow/Backup';
save_folder = '/media/dinesh/M2SSD/Odor Concentration Analysis/Airflow/Sorted';

% %% Get the experiment folder, output folder and backup folder
% expDir = uigetdir(pwd,'Experiment Folder (Trackit save folder)');
% rootDir = fileparts(expDir);
% outDir = uigetdir(rootDir,'Output Folder (Main Data folder)');
% rootDir = fileparts(outDir);
% backupDir = uigetdir(rootDir,'Backup Folder (To move data)');
% rootDir = fileparts(backupDir);
% save_folder = uigetdir(rootDir,'Sorted data folder (sorted data)');

%% Process trackit data
raw_data_matfile = copyTrackitExperimentData(expDir,outDir,backupDir,backup_on,force_rewrite);
all_traj_matfile = extractTrajectories(raw_data_matfile, force_rewrite);
sorted_matfile = extractFirstLandedTrajectories(all_traj_matfile, force_rewrite);

%% Plot treatmentwise trajectory
TrajPlotTreatment(sorted_matfile);
saveSortedTrackItData(sorted_matfile, save_folder);

close all;
