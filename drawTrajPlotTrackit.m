function [] = drawTrajPlotTrackit(pts,cs,varargin)

% Initial Treatment
% If its ov0cm - then mention, else make it 0.
% treatment = 'ov0cm';
rcolor = [1 0 1];
bcolor = [0 0 1];
marker = 'none';
markersize = 5;
markerfill = 'none';

if (isempty(varargin))
    colors = repmat([0 0 1],length(pts)-1,1);        
    decorate = 1;
    tags = 0;
elseif size(varargin,2)==1
    colors = varargin{1};
    decorate = 1;
    tags = 0;
elseif size(varargin,2)==2
    colors = varargin{1};
    decorate = varargin{2};
    tags = 0;  
else 
    colors = varargin{1};
    decorate = varargin{2};
    tags = 1;
    textString =  varargin{3};  
end

if decorate==1
    % Initialize
    plot3(0,0,0);
    hold on;
    grid on;    
end

% Plot trajectory
plot3(pts(:,1),pts(:,2),pts(:,3),'color',colors,'Marker',marker,...
    'MarkerSize',markersize,'MarkerFaceColor',markerfill);
% Mark the start points with filled 'x'
plot3(pts(1,1),pts(1,2),pts(1,3),'color',rcolor,'Marker','x','MarkerFaceColor',rcolor);
% Mark the end points with filled 'o'
plot3(pts(length(pts),1),pts(length(pts),2),pts(length(pts),3),'color',bcolor,'Marker','o','MarkerFaceColor',bcolor);

if (tags==1)
    textString = strjoin(strsplit(textString,'_'),'-');
    text(pts(1,1),pts(1,2),pts(1,3),strcat('   ',char(textString)));
end

% Decorate plot, i.e. add spheres and stuff
if decorate==1
decoratePlot(cs);
end

end

function [] = decoratePlot(cs)
%
%
% Break down treatment
% C = strsplit(treatment,'_');
% treatment = C{1};
% airflow = C{2};

% Initialize sphere
[x,y,z]=sphere;
color=zeros(21,21,3);
row=size(cs,1); 

% Plot other cues
for i = 1:row
    % Plot visual cue sphere
    a=(0.003*x)+cs(i,1);
    b=(0.003*y)+cs(i,2);
    c=(0.003*z)+cs(i,3);
    surf(a,b,c,color);
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