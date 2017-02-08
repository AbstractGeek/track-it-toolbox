function [tags,datenums,varnames,trialnames] = ...
    decipherTrialNamesAndTags (dirList)
% function [] = decipherTrialNamesAndTags (dirList)
% 
% 
% 

%% Defaults
time_count = 6; % Number of str pieces describing time

%% Split strings
splitstrs = cellfun(@strsplit,dirList,repmat({'_'},length(dirList),1),...
    'UniformOutput',false);
tags = cellfun(@(x) strjoin(x(time_count+1:end),'_'),splitstrs,...
        'UniformOutput',false);
datenums = cell2mat(cellfun(@(x) datenum(strjoin(x(1:time_count),'-'),...
    'yyyy-mm-dd-HH-MM-SS'),splitstrs,'UniformOutput',false));
varnames = cellfun(@(x,y) sprintf('e_%s_%s',y,x),tags,...
    cellstr(datestr(datenums,'yyyy_mm_dd')),'UniformOutput',false);
trialnames = cellstr(datestr(datenums,'yyyy_mm_dd_HH_MM_SS'));

end