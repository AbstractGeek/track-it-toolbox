function [] = extractTrajectories(varargin)
% function [] = extractTrajectories()
%
%
%
% Dinesh Natesan, 4th Aug 2014

% Generate plots? for debug mode
generatePlots = 1;
objTag = 'objpositions';

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
rawData = load(fullfile(PathName,FileName));

%% Get object location
names = fieldnames(rawData);
attributes = rawData.(names{~cellfun(@isempty,regexp(names,'attr'))});
objTagLocation = ~cellfun(@isempty,regexp(attributes.tags,objTag));
objPositions = extractObjectLocations(rawData.(attributes.varnames{objTagLocation}));
expname = sprintf('d_%s',datestr(attributes.datenums(1),'yyyy_mm_dd'));

%% Extract trajectories
% Start processing rest of the data
varnames = attributes.varnames;
% Load existing mat file and create necessary directories
if exist(fullfile(PathName,'Trajectories'),'file')==2
    trajData = load(fullfile(PathName,'Trajectories'));
else
    trajData = struct;
end
if ~isdir(fullfile(PathName,'TrajPlots',expname)) && generatePlots == 1
    mkdir(fullfile(PathName,'TrajPlots',expname));
end

% Create Figure
close all;
h = figure();

% Begin appending onto the structure
for j=1:length(varnames)
    % Extract trajectories
    xyzData = extractXYZfromCameraData(rawData.(varnames{j}));
    trajData.(expname).(varnames{j}) = xyzData;
    % Create a figure containing all the trials. Uncomment if necessary
    if (generatePlots == 1)
        colors = colormap(lines(length(xyzData)));                
        
        for k=1:length(xyzData)            
            if k==1
                drawTrajPlotTrackit(xyzData{k}(:,2:4),objPositions,...
                    colors(k,:),1,varnames{j});
            else
                drawTrajPlotTrackit(xyzData{k}(:,2:4),objPositions,...
                colors(k,:),0,varnames{j});
            end
        end
        
        saveas(h,strcat(fullfile(PathName,'TrajPlots',expname,varnames{j}),'.fig'));
        clf(h);
    end
end
close(h);
% Save trajData
save(fullfile(PathName,'Trajectories'),'-struct','trajData');
% Done, exit

end