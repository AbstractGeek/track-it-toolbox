function [trackit_traj_matfile] = extractTrajectories(matfile, force_rewrite)
% function [] = extractTrajectories()
%
%
%
% Dinesh Natesan, 4th Aug 2014

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
        
        if exist(fullfile(currDir, sprintf('%s_cameraData.mat', days{j})),'file') == 2
            trackit_camData = load(fullfile(currDir,...
                sprintf('%s_cameraData.mat', days{j})));
        else
            trackit_camData = struct;
        end
       
        for k=1:length(trials)
            
            outDir = fullfile(currDir,trials{k});
            
            if ~isdir(outDir)
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
            trajFileName = fullfile(outDir, sprintf('%s_%s',days{j},trials{k}));
            TrajPlotTrackitObjects(sortedData, trajFileName);            
            
            % Add entry into log file
            fprintf(logid, '%s: Successfully extracted and plotted %d trajectories\n',...
                sprintf('%s_%s_%s',days{j}(2:end),...
                trackit_data.(treatments{i}).name,trials{k}),length(objects));            
            
            % clear unneccesary vars
            clearvars sortedData cameraData cameraMat;
            
        end
        
        save(fullfile(currDir, sprintf('%s_cameraData.mat', days{j})),...
            '-struct', 'trackit_camData');
        clearvars trackit_camData;
    end
    
    % Save name
    trackit_traj.(treatments{i}).name = trackit_data.(treatments{i}).name;   
    
    % Save trajectories of this treatment as a mat-file
    temp_struct.(treatments{i}) = trackit_traj.(treatments{i}); %#ok<STRNU>
    save(fullfile(rootdir,trackit_data.(treatments{i}).name,...
        sprintf('%s_all_trajectories.mat',trackit_data.(treatments{i}).name)),...
        '-struct', 'temp_struct');       
    clearvars temp_struct;    
end

% Save mat file
save(trackit_traj_matfile,'-struct','trackit_traj');

fprintf(logid, '\nWhole trackit traj data mat file successfully updated and saved\n\n');
fprintf(logid, '\n\nTRAJECTORY EXTRACTION SUCCESSFUL! \n\n');
fclose(logid);

close(h1);
end
