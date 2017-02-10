function [trackit_sorted_data_matfile] = extractFirstLandedTrajectories(matfile, force_rewrite)
% 
% 
% 
% Dinesh Natesan, 25th Aug 2014
% 



if (nargin<1)
    error('extractFirstLandedTrajectories needs a rootdir input to load matfiles');
elseif nargin == 1
    force_rewrite = 0;
end

% Out file
trackit_sorted_data_mat = 'trackit_sorted_data.mat';

rootdir = fileparts(matfile);
trackit_traj = load(matfile);

trackit_sorted_data_matfile = fullfile(rootdir, trackit_sorted_data_mat);

if exist(trackit_sorted_data_matfile,'file') == 7
    trackit_sorted_data = load(trackit_sorted_data_matfile);   
else 
    trackit_sorted_data = struct;
end

%% Extract First Landed Trajectory

treatments = fieldnames(trackit_traj);

for i=1:length(treatments)
    
    days = fieldnames(trackit_traj.(treatments{i}));
    days(ismember(days, {'name'})) = [];
    
    csvfolder = fullfile(rootdir,trackit_traj.(treatments{i}).name);
    trajfolder = fullfile(rootdir,trackit_traj.(treatments{i}).name,'Traj-plots (Daywise)');
    
    if ~isdir(trajfolder)
        mkdir(trajfolder);
    end
    
    for j=1:length(days)
        
        currDir = fullfile(rootdir,trackit_traj.(treatments{i}).name,...
                'Sorted-Data',days{j});  
        
        % Check if log_file exists
        log_file = fullfile(currDir, 'trackit_first_landed_extraction.log');
        
        if (exist(log_file, 'file')==7) && (force_rewrite==0)
           % Skip extraction (already done)
           continue;
        end
        
        % Append events to log file
        logid = fopen(log_file, 'a');
        fprintf(logid,'\n#################### %s ####################\n\n',...
            datestr(now));
        
        % Perform object extraction
        [cs,trials] = ...
            extractObjectLocations(trackit_traj.(treatments{i}).(days{j}));
        
        if isempty(trials) && length(cs) == 1
            fprintf('Found %d object trials in day %s of treatment %s. Unfamiliar experiment protocol. Skipping day\n',...
                cs, days{j}(3:end), treatments{i});
            fprintf(logid, 'Found %d object trials in day %s of treatment %s. Unfamiliar experiment protocol. Skipping day\n',...
                cs, days{j}(3:end), treatments{i});
            fclose(logid);
            continue; 
        elseif isempty(trials) && (length(cs) == 2)
            fprintf('Found %d objects in day %s of treatment %s. Unfamiliar number of objects. Skipping day\n',...
                cs(2), days{j}(3:end), treatments{i});
            fprintf(logid, 'Found %d objects in day %s of treatment %s. Unfamiliar number of objects. Skipping day\n',...
                cs(2), days{j}(3:end), treatments{i});
            fclose(logid);
            continue;           
        end
        
        % Add to log
        fprintf(logid,'Object Extraction successful!. Beginning first landed trajectory extraction\n\n');        

            for k=1:length(trials)
                
                [trackit_sorted_data.(treatments{i}).(days{j}).(trials{k}), selected_fly] = ...
                    getFirstLanded(trackit_traj.(treatments{i}).(days{j}).(trials{k}), cs);                
                
                % Save data as csv and copy to main folder
                writetable(trackit_sorted_data.(treatments{i}).(days{j}).(trials{k}),...
                    fullfile(currDir, sprintf('%s_%s_%s.csv',days{j}(3:end),...
                    trackit_traj.(treatments{i}).name,trials{k})));
                copyfile(fullfile(currDir, sprintf('%s_%s_%s.csv',days{j}(3:end),...
                    trackit_traj.(treatments{i}).name,trials{k})), csvfolder);
                
                % Add to log
                fprintf(logid, '%s: %s selected as the first landed fly.\nTrajectory successfully saved as csv\n\n',...
                    sprintf('%s_%s_%s',days{j}(3:end),...
                trackit_traj.(treatments{i}).name,trials{k}), selected_fly);
            
            end
            
            % Plot trajectories
            TrajPlotTrackitObjects(trackit_sorted_data.(treatments{i}).(days{j}),...
                fullfile(trajfolder, sprintf('%s_%s.fig',days{j}(3:end),...
                trackit_traj.(treatments{i}).name)));            
                
            fprintf(logid, 'FIRST LANDED TRAJECTORY EXTRACTION SUCCESSFUL\n\n');        
            fclose(logid);        
    end
    
    % Save name
    trackit_sorted_data.(treatments{i}).name = trackit_traj.(treatments{i}).name;   
    
    % Save first landed trajectories of this treatment as a mat-file
    temp_struct.(treatments{i}) = trackit_sorted_data.(treatments{i}); %#ok<STRNU>
    save(fullfile(rootdir,trackit_traj.(treatments{i}).name,...
        sprintf('%s_sorted_data.mat',trackit_traj.(treatments{i}).name)),...
        '-struct', 'temp_struct');
    clearvars temp_struct;    
        
end

% Save mat file
save(trackit_sorted_data_matfile,'-struct','trackit_sorted_data');

end


function [firstLanded, selected_fly] = getFirstLanded(trialData, cs)
% [firstLanded, selected_obj] = getFirstLanded(trialData, cs)
% 
% Inputs: Trial Data: Trial Data structure with multiple trajectors
%         cs: Object center (X,Y,Z)
% 
% Outputs: firstLanded: Table with trajectory of the first landed fly along
%                       with the object centers used for selection.
%          selected_obj: The selected fly based on first landed criterion.
% 
% Dinesh Natesan
% 10th Feb, 2017

objects = fieldnames(trialData);
minimum_dist = nan(size(objects));
land_time = nan(size(objects));
start_dist = false(size(objects));

dist_threshold = 0.08;      % 8 cms
dist_roundoff = 2;          % 1 cm round off is good enough

for i = 1:length(objects)
    minimum_dist(i) = floorn(sqrt(sum([trialData.(objects{i}).X(end) - cs(1),...
        trialData.(objects{i}).Y(end) - cs(2),...
        trialData.(objects{i}).Z(end) - cs(3)].^2)),dist_roundoff);
    land_time(i) = trialData.(objects{i}).time(end);
    start_dist(i) = (sqrt(sum([trialData.(objects{i}).X(1) - cs(1),...
        trialData.(objects{i}).Y(1) - cs(2),...
        trialData.(objects{i}).Z(1) - cs(3)].^2))) >= dist_threshold;
end

obj_table = sortrows(table(objects(start_dist), minimum_dist(start_dist),...
    land_time(start_dist),...
    'VariableNames',{'objects','minimum_dist','land_time'}),...
    {'minimum_dist'},{'ascend'});
selected_fly = obj_table.objects{1};

% Create table and reset time
firstLanded =  trialData.(selected_fly);
firstLanded.time = firstLanded.time - firstLanded.time(1);
table_rows = length(firstLanded.time);

% Obtain cs headers (handles multiple cs)
xyz_tag = repmat({'_X', '_Y', '_Z'},size(cs,1),1);
cs_strings = repmat(arrayfun(@(x) sprintf('cs%02d',x),...
    (1:size(cs,1)),'UniformOutput', false)',1,3);
cs_headers = reshape(strcat(reshape(cs_strings,3*size(cs,1),1),...
    reshape(xyz_tag,3*size(cs,1),1)),size(cs,1),3);

% reshape and format cs and its headers
reshaped_cs_headers = reshape(cs_headers',size(cs_headers,1)*size(cs_headers,2),1)';
reshaped_cs = reshape(cs',size(cs,1)*size(cs,2),1)';

% Create a temp table with cs
temp_table = array2table(nan(table_rows,3*size(cs,1)),'VariableNames',reshaped_cs_headers);
temp_table{1,:} =  reshaped_cs;

% add it to the main table
firstLanded = [firstLanded, temp_table];

end


function [Xr] = floorn(X,n)
% function [X] = floorn(X,n)
% 
% Floors the input at the nth decimal place (similar to round at nth
% decimal point)
% 
% Dinesh Natesan
% 10th Feb 2015

Xr = floor(X*10^n)*10^-n;
end