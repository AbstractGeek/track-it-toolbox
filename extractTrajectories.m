function [trackit_traj_matfile] = extractTrajectories(matfile, force_rewrite)
% function [trackit_traj_matfile] = extractTrajectories(matfile, force_rewrite)
%
% Sift through all the trajectories captured by the trackit system and
% extract them. Save individual trajectories as csv files and combined all
% the data into a matfile.
% 
% Inputs:
%   matfile: 
%       File location of the matfile containing raw data from trackit 
%       experiments(output of copyTrackitExperimentData function).
%   force_rewrite: 
%       Flag to rewrite the data even if the analysis has been done before.
% 
% Outputs:
%   trackit_traj_matfile: 
%       a mat file containing all the trajectories from all trackit
%       experiments (sorted based on treatments and experiment days). 
% 
% Dinesh Natesan 

if (nargin<1)
    error('extractTrajectories needs a matfile input to load raw data');
elseif nargin == 1
    force_rewrite = 0;
end

% Out file
trackit_traj_mat = 'trackit_all_trajectories.mat';

rootdir = fileparts(matfile);
trackit_data = load(matfile);

trackit_traj_matfile = fullfile(rootdir, trackit_traj_mat);

if exist(trackit_traj_matfile,'file') == 2
    trackit_traj = load(trackit_traj_matfile);   
else 
    trackit_traj = struct;
end

% Append events to a log file
log_file = fullfile(rootdir, 'trackit_all_traj_extraction.log');
logid = fopen(log_file, 'a');
fprintf(logid,'\n#################### %s ####################\n\n',datestr(now));

% Create figure
h1 = figure();

%% Extract trajectories
treatments = fieldnames(trackit_data);

for i=1:length(treatments)
    
    days = fieldnames(trackit_data.(treatments{i}));
    days(ismember(days, {'name'})) = [];
    
    for j=1:length(days)
        
        trials = fieldnames(trackit_data.(treatments{i}).(days{j}));
        currDir = fullfile(rootdir,trackit_data.(treatments{i}).name,...
                'Sorted-Data',days{j});             
        figDir = fullfile(rootdir,trackit_data.(treatments{i}).name,...
            'Traj-plots',days{j});
        
        if ~isfolder(figDir)
            mkdir(figDir);
        end
        
        if exist(fullfile(currDir, sprintf('%s_cameraData.mat', days{j})),'file') == 2
            trackit_camData = load(fullfile(currDir,...
                sprintf('%s_cameraData.mat', days{j})));
        else
            trackit_camData = struct;
        end
       
        for k=1:length(trials)
            
            outDir = fullfile(currDir,trials{k});
            
            if ~isfolder(outDir)
                % Folder doesn't exist, hence data not processed
                mkdir(outDir);
            elseif (force_rewrite == 0)
                % Folder exists (data processed once) and no rewrite.
                continue;                
            end
                        
            % Extract objects
            [sortedData, cameraData, cameraMat] = ...
                extractXYZfromCameraData(trackit_data.(treatments{i}).(days{j}).(trials{k}));
            trackit_traj.(treatments{i}).(days{j}).(trials{k}) = sortedData;
            trackit_camData.(days{j}).(trials{k}) = cameraData;
            
            % Save cameraMat as csvs
            csvheaders = cameraMat.headers;
            objects = fieldnames(cameraData);            
            for l = 1:length(objects)
                writeCSV(fullfile(outDir,sprintf('%s.csv',objects{l})),...
                    csvheaders, cameraMat.(objects{l}));
            end
            
            % Create trajectories of each object
            trajFileName = fullfile(figDir, sprintf('%s_raw',trials{k}));
            TrajPlotTrackitObjects(sortedData, trajFileName);            
            
            % Add entry into log file
            fprintf(logid, '%s: Successfully extracted and plotted %d trajectories\n',...
                sprintf('%s_%s_%s',days{j}(2:end),...
                trackit_data.(treatments{i}).name,trials{k}),length(objects));            
            
            % clear unneccesary vars
            clearvars sortedData cameraData cameraMat;
            
        end
        
        save(fullfile(currDir, sprintf('%s_cameraData.mat', days{j})),...
            '-struct', 'trackit_camData','-v7.3');
        clearvars trackit_camData;
    end
    
    % Save name
    trackit_traj.(treatments{i}).name = trackit_data.(treatments{i}).name;   
    
    % Save trajectories of this treatment as a mat-file
    temp_struct.(treatments{i}) = trackit_traj.(treatments{i}); %#ok<STRNU>
    save(fullfile(rootdir,trackit_data.(treatments{i}).name,...
        sprintf('%s_all_trajectories.mat',trackit_data.(treatments{i}).name)),...
        '-struct', 'temp_struct','-v7.3');       
    clearvars temp_struct;    
end

% Save mat file
save(trackit_traj_matfile,'-struct','trackit_traj','-v7.3');

fprintf(logid, '\nWhole trackit traj data mat file successfully updated and saved\n\n');
fprintf(logid, '\n\nTRAJECTORY EXTRACTION SUCCESSFUL! \n\n');
fclose(logid);

close(h1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = TrajPlotTrackitObjects(sortedData, filename)
% function [] = TrajPlotTrackitObjects(sortedData, filename)
%
% A simple function to plot all the objects detected by the trackit system
% in one experiment trial.
% 
% Inputs: 
%   sortedData: A structure with tables of object trajectories.
%   filename: Filename of the plotted 3D figure.
% 
% Dinesh Natesan 

objects = fieldnames(sortedData);

% Initialize
plot3(0,0,0,'HandleVisibility','off');
hold on;
grid on;

% Get colors
colors = lines(length(objects));

% Plot all data
for i = 1:length(objects)    
    plot3(sortedData.(objects{i}).X, sortedData.(objects{i}).Y, sortedData.(objects{i}).Z,...
        'Color', colors(i,:), 'DisplayName', objects{i})    
end

legend show;

% Save data
savefig(filename);

% Clear figure
clf(gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sortedData, cameraData, cameraMat] = extractXYZfromCameraData(rawData)
% function [sortedData, cameraData, cameraMat] = extractXYZfromCameraData(rawData)
% Extracts out time, X, Y, Z from the input data for all cameras (camera
% number independent) and averages X, Y, Z of cameras in the same time point
% for each object in the dataset.
% 
% Input: rawData - table generated from Filtered file from trackit
% 
% Output: 
% sortedData - contains a table for each object tracked by trackit.
% table contains time and averaged X, Y, Z from all camera views for each
% time point.
% 
% cameraData - contains a table for each object tracked by trackit.
% table containing time, X, Y, Z for each camera. This is the raw data used
% to generate the table in sortedData
% 
% cameraMat - same data as the cameraData table except save in mat files so
% as to allow easy conversion to csv files.
% 
% Dinesh Natesan 
% 9th Feb 2017

% Defaults
t_round_dec_pt = 4;       % max_frame_rate = 120; 0.0001 resolution should be more than enough

% Sort by objects
obj = sort(unique(rawData.object));
cameras = sort(unique(rawData.camera));
sortedData = struct;
cameraData = struct;
cameraMat = struct;

for i = 1:length(obj)    
    % Extract object indexes and associated time values
    objInd = (rawData.object== obj(i));
    time = sort(unique(round(rawData.time(objInd),t_round_dec_pt)));
    
    % Initialize camera variables
    camera_time = nan(length(time),length(cameras));
    camera_X = nan(length(time),length(cameras));
    camera_Y = nan(length(time),length(cameras));
    camera_Z = nan(length(time),length(cameras));
    
    % Extract camera specific values
    for j=1:length(cameras)
        % Extract camera indices of this object
        camera_ind = (rawData.camera == cameras(j)) & objInd;
        % Extract time, X, Y, Z
        temp_time = round(rawData.time(camera_ind),t_round_dec_pt);
        temp_X = rawData.X(camera_ind);
        temp_Y = rawData.Y(camera_ind);
        temp_Z = rawData.Z(camera_ind);
        % Find intersection and copy data
        [~,ia,~] = intersect(time,temp_time);
        camera_time(ia,j) = temp_time;
        camera_X(ia,j) = temp_X;
        camera_Y(ia,j) = temp_Y;
        camera_Z(ia,j) = temp_Z;
        % Clear temp variables
        clearvars temp_time temp_X temp_Y temp_Z;
    end
    
    % Sanity check
    if ~all(time == nanmean(camera_time,2))
       error('Sanity Check didnt work in extractXYZfromCameraData');
    end
    
    % Average and fill sorted data table    
    sortedData.(sprintf('obj%04d',obj(i))) = ...
        table(time,nanmean(camera_X,2),nanmean(camera_Y,2),nanmean(camera_Z,2),...
        'VariableNames',{'time','X','Y','Z'});  
    
    % Generate struct of sorted data before averaging
    cameraData.(sprintf('obj%04d',obj(i))) = table(camera_time, camera_X, camera_Y, camera_Z);    
    cameraMat.(sprintf('obj%04d',obj(i))) = [camera_time, camera_X, camera_Y, camera_Z];
   
end

% Create headers
camera_names = arrayfun(@(x) sprintf('cam%01d', x), cameras, 'UniformOutput', false);
cameraMat.headers = [strcat(camera_names,'_time')',strcat(camera_names,'_X')',...
    strcat(camera_names,'_Y')',strcat(camera_names,'_Z')'];

end
