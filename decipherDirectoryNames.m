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