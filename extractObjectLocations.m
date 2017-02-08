function [cs] = extractObjectLocations(cameradata)
% function [objects] = extractObjectLocations(csvfile)
% 
% 
% 

%% Get object xyz from cameradata
pts = extractXYZfromCameraData(cameradata);        

%% Extract the object center from the location
cs = NaN(length(pts),3);
for i=1:length(pts)
    curr = pts{i};
    cs(i,:) = mean(curr(:,2:4));
end


end