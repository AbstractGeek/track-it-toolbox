function [xyzData] = extractXYZfromCameraData(camData)
% function [xyzData] = extractTrajectories(camData)
% Pulls out the Time,X,Y,Z from the inputted data for both cameras
% Averages out the X,Y,Z for both cameras and outputs the result
% Ignores the single frames
% Output format - each cell has one trajectory
%               - Contents of a cell [Time,X,Y,Z]
% 
% 
% Dinesh Natesan, 4th Aug 2014
% 

% Sort by objects
obj = unique(camData(:,2));
xyzData = cell(length(obj),1);
delList = [];
for i = 1:length(obj)
    objInd = (camData(:,2)== obj(i));
    tempData = [camData(objInd,1),camData(objInd,3),camData(objInd,16:18)];
    toMerge=[false;diff(tempData(:,2))<=1e-4]; %mergeData that have been measured within a 100us interval
    % Ignoring single data and picking ones that were detected by both
    % cameras
    mergedTime=(tempData(toMerge,2)+tempData([toMerge(2:end);false],2))/2;
    mergedX=(tempData(toMerge,3)+tempData([toMerge(2:end);false],3))/2;
    mergedY=(tempData(toMerge,4)+tempData([toMerge(2:end);false],4))/2;
    mergedZ=(tempData(toMerge,5)+tempData([toMerge(2:end);false],5))/2;
    xyzData{i,1} = [mergedTime,mergedX,mergedY,mergedZ];
    if isempty(xyzData{i,1})
        delList = [delList;i];
    end
    % Done - leave
end
% Delete all empty cell arrays
xyzData(delList)=[];
% Done, leave
end