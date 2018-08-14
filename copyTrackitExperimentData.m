function [raw_data_matfile] = copyTrackitExperimentData(expDir,outDir,backupDir,backup_on,force_rewrite)
% function [raw_data_matfile] = copyTrackitExperimentData(expDir,outDir,backupDir,force_rewrite)
%
% Copy all experiment data from experiment directory to output directory,
% while backing the key files in the backup directory. The code
% automatically identifies the date, the treatment and the tag for each
% experiment trial and uses it for sorting. Apart from copying and backing
% up experiment data, the code also saves the data as a structure sorted
% based on treatment and experiment day. It also automatically links the
% trials to the object positions and adds it appropriately (see note
% below). The raw data is finally saved as a mat_file and its location is
% outputted.
%
% Note: Folder names of experiment trials show have this format:
%       "YYYY_MM_DD_treatmentname_tag".
%       tag can be either trial number (automatically assigned if not
%       present) or object positions which are required for analysis. T
%       The code requires one object position tag for every treatment
%       perform per day in order to consider it a valid experiment.
%
% Inputs:
%       expDir: Directory containing folders of experiments.
%       outDir: Directory to copy (and sort) experiment data.
%       backupDir: Directory to backup the experiment data.
%       force_rewrite:  Flag to rewrite the data even if the
%                       analysis has been done before.
%
% Outputs:
%   raw_data_matfile:
%       a mat file containing the raw data from experiments (sorted based
%       on treatments and experiment days).
%
% Dinesh Natesan


minTrajLength = 30;      % frames
trackit_data_mat = 'trackit_raw_data.mat';

if (nargin<3)
    error('extractTrajectories needs a matfile input to load raw data');
elseif nargin == 3
    backup_on = 0;
    force_rewrite = 0;
end

%% Let's rock and roll?
% Obtain a list of subfolders
listOfContents = dir(expDir);
isubDir = [listOfContents(:).isdir];
dirList = {listOfContents(isubDir).name}';
dirList(ismember(dirList,{'.','..'})) = [];

% Append events to a log file in the output directory
log_file = fullfile(outDir, 'trackit_data_copy.log');
logid = fopen(log_file, 'a');
logcleanup = @() fclose(logid);
fprintf(logid,'\n#################### %s ####################\n\n',datestr(now));

% Clean up directory names
fprintf(logid, 'Checking directory names\n');
dirTable = decipherDirectoryNames (dirList, logid);
dirNum = size(dirTable,1);

% Save raw data mat file location
raw_data_matfile = fullfile(outDir,trackit_data_mat);

% Load mat file directory
if (exist(raw_data_matfile,'file') == 2)
    trackit_data = load(raw_data_matfile);
else
    trackit_data = struct;
end

%% Begin by copying the headers
% open a waitbar
h = waitbar((1-1)/(dirNum-1),strcat('Processing: Trial ',num2str(1),'-',num2str(dirNum),'...'));
pos = get(h,'Position');
set(h,'Position',[pos(1) pos(2)+2.*pos(4) pos(3) pos(4)]);

fprintf(logid, '\nBeginning copy:\n');

% Begin loading data into a matfile.
% Copy the csv file as a backup
for i=1:dirNum
    waitbar((i-1)/(dirNum-1),h,strcat('Processing: Trial ',num2str(i),'-',num2str(dirNum),'...'));
    currDir = dirTable.expDirs{i};
    currFile = fullfile(expDir,currDir,dirTable.inFile{i});
    if (exist(currFile, 'file') ~= 2)
        % file doesn't exist
        % show a warning and continue
        fprintf('File: %s doesn''t exit. Skipping\n', currDir);
        fprintf(logid, 'File: %s doesn''t exist. Skipping\n', currDir);
        continue;
    end
    
    % Check if current file already exists in the data and the backup
    % folder
    if ((exist(fullfile(fullfile(outDir, dirTable.treatment{i}, 'Raw-Data'),...
            sprintf('Filtered_%s.csv',currDir)), 'file') == 2) && ...
            (force_rewrite == 0))
        % Folder has already been processed.
        continue;
    end
    
    
    currData = readtable(currFile);
    if ~isempty(currData) && (size(currData,1) > minTrajLength)
        % Save the trajectory data into a mat folder
        treatment_name = dirTable.treatment{i};
        treatment_name(ismember(treatment_name, '(.)')) = []; % make it valid variable name
        trackit_data.(treatment_name).(sprintf('d_%s',datestr(dirTable.dateNums(i),'yyyy_mm_dd'))).(dirTable.trialDetails{i}) = currData;
        
        % Add a treatment complete name string if it doesn't exist
        if ~isfield(trackit_data.(treatment_name), 'name')
            trackit_data.(treatment_name).name = dirTable.treatment{i};
        end
        
        % Copy filtered data into the treatment folder
        treatmentDir = fullfile(outDir, dirTable.treatment{i}, 'Raw-Data');
        if ~isdir(treatmentDir)
            % Add a entry into log file
            mkdir(treatmentDir);
        end
        copyfile(currFile,fullfile(treatmentDir,sprintf('Filtered_%s.csv',currDir)));
        
        if backup_on
            % Move the whole folder to the backup folder
            treatmentDir = fullfile(backupDir, dirTable.treatment{i});
            if ~isdir(treatmentDir)
                % Add a entry into log file
                mkdir(treatmentDir);
            end
            copyfile(fullfile(expDir,currDir),fullfile(treatmentDir,currDir));
        end
        
        fprintf(logid, '%s: Successfully processed\n', currDir);
    else
        fprintf('%s has no data: Not processed\n', currDir);
        fprintf(logid, '%s has no data: Not processed\n', currDir);
    end
    
end
close(h);

% Save mat file
save(raw_data_matfile,'-struct','trackit_data','-v7.3');

fprintf(logid, '\nWhole trackit_data mat file successfully updated and saved\n\n');

fprintf(logid, 'Updating individual treatment mat files:\n');
% Save individual treatments
treatments = fieldnames(trackit_data);
for i = 1:length(treatments)
    temp_struct.(treatments{i}) = trackit_data.(treatments{i}); %#ok<STRNU>
    save(fullfile(outDir,trackit_data.(treatments{i}).name,sprintf('%s_raw_data.mat',trackit_data.(treatments{i}).name)),...
        '-struct', 'temp_struct');
    fprintf(logid, '%s: Successful \n', trackit_data.(treatments{i}).name);
    clearvars temp_struct
end

fprintf(logid, '\n\nCOPY SUCCESSFUL! \n\n');

end

%%% Subfunctions necessary for proper running %%%

function [dirTable] = decipherDirectoryNames (dirList, logid)
% function [dirTable] = decipherDirectoryNames (dirList)
%
% The input dirList is a set of directories in the Trackit Data folder.
% The output dirTable is a table with all the important attributes
% extracted from the directory name.
%
% The folders should have the following naming convention:
% 'yyyy_mm_dd_HH_MM_SS_$treatmentname_$trialdetail'
%
% $treatmentname should be an abbreviation for the experiment treatment. It
% can contain most of the alphanumeric characters that is generally allowed
% in a folder name except '_'
%
% $trialdetail should contain key details about the experimental trial
% (example trial number). It too can contain most of the alphanumeric
% characters that is generally allowed in a folder name except '_'
%
% Dinesh Natesan
% 8th Feb 2017


% defaults
time_count = 6; % Number of str pieces describing time

% split strings
split_names = cellfun(@strsplit,dirList,repmat({'_'},length(dirList),1),...
    'UniformOutput',false);
split_lengths = cellfun(@length, split_names);

% remove invalid directories and print them
valid_dir_names = split_lengths > 6;
expDirs = dirList(valid_dir_names);
if (any(not(valid_dir_names)))
    fprintf('Failed import folder %s. Naming does not match usual convention \n', dirList{not(valid_dir_names)});
    fprintf(logid, 'Failed import folder %s. Naming does not match usual convention \n', dirList{not(valid_dir_names)});
end

% extract datenums and treatment names (required attributes)
dateNums = cell2mat(cellfun(@(x) datenum(strjoin(x(1:time_count),'-'),'yyyy-mm-dd-HH-MM-SS'),...
    split_names(valid_dir_names),'UniformOutput',false));
treatment = cellfun(@(x) x{7}, split_names(valid_dir_names),'UniformOutput',false);
trialDetails = cell(size(treatment));

% extract trial details and give it default value if absent
trialDetails(split_lengths>=8) = cellfun(@(x) x{8}, split_names(split_lengths>=8),'UniformOutput',false);
trialDetails(cellfun(@isempty,trialDetails)) = arrayfun(@(x) sprintf('trial%03d',x),...
    find(split_lengths<8 & split_lengths>6),'UniformOutput',false);

% Generate input and output files
inFile = arrayfun(@(x) sprintf('Filtered_%s.csv',datestr(x,'yyyy_mm_dd_HH_MM_SS')),...
    dateNums,'UniformOutput',false);
outFile = arrayfun(@(x,y) sprintf('%s_%s.csv',datestr(x,'yyyy-mm-dd'),y{1}),...
    dateNums,trialDetails,'UniformOutput',false);

% Create table
dirTable = table(expDirs, dateNums, treatment, trialDetails, inFile, outFile);

end