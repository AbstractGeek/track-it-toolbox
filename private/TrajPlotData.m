function [] = TrajPlotData(cs, treatment, sortedData, filename)
% function [] = TrajPlotData(cs, treatment, sortedData, filename)
% 
% Plot the trajectory data for the sortedData along with the centers of the
% object (based on treatment tag). 
% 
% Inputs:
%   cs: centers of objects in the treatment.
%   treatment: treatment tag of the experiment.
%   sortedData: A structure containing tables of trajectories.
%   filename: output figure filename
% 
% Output: 
%   Saved 3d trajectory figure.
% 
% Dinesh Natesan 
% Last updated: 12th Oct 2017

objects = fieldnames(sortedData);

% Initialize
plot3(0,0,0,'HandleVisibility','off');
hold on;
grid on;

% Get colors
colors = lines(length(objects));

% Plot all data
for i = 1:length(objects)    
    plot3(sortedData.(objects{i}).X, sortedData.(objects{i}).Y, sortedData.(objects{i}).Z,...
        'Color', colors(i,:), 'DisplayName', objects{i})    
end

legend show;

% Decorate plot
decoratePlot(cs, treatment)

% Save data
savefig(filename);

% Clear figure
clf(gcf);
end

function [] = decoratePlot(cs, treatment)
%
%

% Defaults
sphere_radius = 0.003;
cylinder_radius = 0.00025;
cylinder_height = 0.005;
redcolor = [0.8980 0 0];
palepinkcolor = [1.0000 0.8118 0.8627];

% Initialize shape
[xs,ys,zs]=sphere;
[xc,yc,zc]=cylinder(cylinder_radius);

% Initialize color
black_sphere = zeros(21,21,3);
red_sphere = reshape(repmat(redcolor,21*21,1),21,21,3);
coral_sphere = reshape(repmat(palepinkcolor,21*21,1),21,21,3);
black_cyl= zeros(2,21,3);
red_cyl = reshape(repmat(redcolor,2*21,1),2,21,3);
palepink_cyl = reshape(repmat(palepinkcolor,2*21,1),2,21,3);

% Assign table for object shape and color
obj_det = table(cell(size(cs,1),1),cell(size(cs,1),1),cell(size(cs,1),1),...
    'VariableNames',{'Shape', 'Color', 'EdgeColor'});

% Cleanup treatment name
treatment = treatment(isletter(treatment)); 
if rem(length(treatment),2) ~= 0
    % Unknown experimental treament
    % Set axis and return
    v=axis;
    pmax=round(max(v));
    pmin=round(min(v));
    % Make axis equal and square
    v1=v; % for now - add a smarter critera later
    axis on;
    axis manual;
    axis equal;
    axis(v1);
    axis vis3d;
    return;
end
obj_num = length(treatment)/2;

% Assign shape based on treatment
low_contrast = treatment(1:obj_num)=='c';
high_contrast = treatment(1:obj_num)=='V';
obj_det.Shape(low_contrast) = {[xc,yc,-zc.*cylinder_height]};    % Low contrast cylinder
obj_det.Shape(high_contrast) = {[xs,ys,zs].*sphere_radius};   % High contrast sphere

% Obtain color characteristics of the treatment
no_odor = ismember((1:obj_num),(find(treatment=='n')-obj_num));
low_odor = ismember((1:obj_num),(find(treatment=='L')-obj_num));
high_odor = ismember((1:obj_num),(find(treatment=='H')-obj_num));

% Look for all combinations of vision and odor
% vision + no odor
obj_det.Color(low_contrast & no_odor) = {black_cyl};
obj_det.Color(high_contrast & no_odor) = {black_sphere};
obj_det.EdgeColor(no_odor) = {[0,0,0]};
% vision + low odor
obj_det.Color(low_contrast & low_odor) = {palepink_cyl};
obj_det.Color(high_contrast & low_odor) = {coral_sphere};
obj_det.EdgeColor(low_odor) = {palepinkcolor};
% vision + high odor
obj_det.Color(low_contrast & high_odor) = {red_cyl};
obj_det.Color(high_contrast & high_odor) = {red_sphere};
obj_det.EdgeColor(high_odor) = {redcolor};

% Plot objects
for i = 1:size(cs,1)
    shape = obj_det.Shape{i};
    a= shape(:,1:21)+cs(i,1);
    b= shape(:,22:42)+cs(i,2);
    c= shape(:,43:63)+cs(i,3);
    surf(a,b,c,obj_det.Color{i},'EdgeColor',obj_det.EdgeColor{i});
end
    
% Set axis
v=axis;
pmax=round(max(v));
pmin=round(min(v));
% Make axis equal and square
v1=v; % for now - add a smarter critera later
axis on;
axis manual;
axis equal;
axis(v1);
axis vis3d;

end
