%% Project: 
% Date: 29/01/22

%% CALCULATION AND REPRESENTATION OF BEZIER CURVES
% Using two diferent methods, visualizing first the drawing of a curve
% using De Casteljau's algorithm, and then the simple calculation of the
% curve.

clear
close all

%% General setup 
% To start: define a points matrix "P" with first row representing x
% coordinate, and second row representing y coordinate.
% Use sample cases below, or define own vector such as:
%   P = [ x1   x2  ...   xn ;
%         y1   y2  ...   yn ];
% Or else activate the point selectino polygon

% Define number of steps desired
steps = 30;

% For animations, set to 1 (will save gif in folder
animation = 0;
imsize = 400;

% Define margins for graph
margin = 0.1;


%% SAMPLE CASES (all points between 0 and 1)
% select 1,2, or 3 for first, second, or third order curves (respectively)
% select 4,5,6 or 7 for higher order (more complex) curves
% select case 8 for a curve with cusp
P = sample_cases(5);

%% POINT SELECTION POLYGON

Fdecas = figure(1); x0=200; y0=200; width=imsize; height=0.93*imsize;
set(1,'position',[x0,y0,width,height])
axis([0 1 0 1]);
title('Click to place points, press "Enter" when finished')
[x, y] = getpts();
clf(1,'reset')
P = [x'; y']';

%% Point checker
if size(P,2)<2
    error('At least two points must be defined')
elseif size(P,1)>size(P,2)
    P = P';
end

%% Preliminary definitions

% Order of curve
n = length(P)-1;

% Time vector
step = 1/steps;
tvec = linspace(0,1,steps);

%% DECASTELJAU ALGORITM ANIMATION

Fdecas = figure(1); x0=200; y0=200; width=imsize; height=0.93*imsize;
set(Fdecas,'position',[x0,y0,width,height])
hold on
title(sprintf('Bezier curve of order %i with De Casteljau algorithm', n));
init(P,margin);

Bdecas = B_casteljau(P,tvec);
im = cell(1,steps);

for i=1:steps
    t = tvec(i);
    
    % Delete the plotted temporary lines and dots
    if i>1
        delete(linemarker)
        delete(dotmarker)
    end
    
    % Plot the Bezier curve from 0 to the calculated instant
    plot(Bdecas(1,(1:i)),Bdecas(2,(1:i)),'k','LineWidth',3);
    
    % Plot the intermediate points
    [linemarker,dotmarker] = intermediate(P,t);
    
    if animation==1
        frame = getframe(Fdecas);
        im{i} = frame2im(frame);
    end
end

if animation==1
    writegif('decas.gif',im,steps,2/steps);
end

%% CLASSICAL CALCULATION OF BEZIER CURVE
% This implementation is much quicker, as the recursive nature adds one
% loop per point, as opposed to one calculatio for each time step and at
% each point.

% Call function
Bclas = B_classic(P,tvec);

% Plot results
Fclas = figure(2); x0=200+imsize; y0=200; width=imsize; height=0.93*imsize; hold on
set(Fclas,'position',[x0,y0,width,height])
title(sprintf('Bezier curve of order %i with classic implementation', n));
init(P,margin);
plot(Bclas(1,:),Bclas(2,:),'.b','LineWidth',3);


%% BEZIER CURVES OF BEZIER CURVES

% After a Bezier curve is calculates, use the points of that curve to
% calculate a second-iteratio bezier curve, and so on

% repeat = 10*n;
% Brep = B_classic(P,tvec);
% 
% figure(3); hold on
% title(sprintf('Bezier curve of Bezier curve, %i times', repeat));
% init(P,margin);
% c = parula(repeat);
% 
% for i=1:repeat
%     Brep = B_classic(Brep,tvec);
%     plot(Brep(1,:),Brep(2,:),'k','LineWidth',3,'Color',c(i,:));
%     
%     if animation==1
%         frame = getframe(Fdecas);
%     end
% end
