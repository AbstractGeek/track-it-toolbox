function [] = extractFirstLandedTrajectories(varargin)
% 
% 
% 
% Dinesh Natesan, 25th Aug 2014
% 

if isempty(varargin)
    % Input mat-file and the load it.
    [FileName,PathName,~] = uigetfile('*.mat');
else
    % Format filename
    [PathName,FileName,ext] = fileparts(varargin{1});
    FileName = strcat(FileName,ext);
    clear ext;
end
% Load Data
TrajData = load(fullfile(PathName,FileName));
% ObjectData = load(fullfile(PathName,''));
% FieldNames
names = fieldnames(TrajData);




end
