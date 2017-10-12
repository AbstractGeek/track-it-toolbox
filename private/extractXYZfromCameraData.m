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