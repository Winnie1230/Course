clear all; close all; clc;

%% create map

axis([-50 50 -15 250])
axis equal
% pos = [10 5 ; -45 30; 5 45 ; -25 60];
% len = [20 20; 30 10; 18 12 ; 15 20];
pos = [5 40 ; -30 90; 5 140 ; -30 190];
len = [25 15; 25 15; 25 15 ; 25 15];

rec_center(:,1) = pos(:,1)+len(:,1)/2;
rec_center(:,2) = pos(:,2)+len(:,2)/2;
Q = ((len(:,1).^2 + len(:,2).^2).^(1/2)/2)';
% Q = Q*5;
% Q = max(len,[],2)
% max_Q = 10;
% min_Q = 5;
% Q = min_Q + (max_Q-min_Q)/(max(Q)-min(Q)).*(Q-min(Q))

for i = 1:length(pos)
    rectangle('Position',[pos(i,:) len(i,:)])
end
hold on
% iter = 200;
% TQ = -500;
% target = [0;90];
% TrajComputation(rec_center,Q,target,TQ,iter,pos,len)

%% car model
center = [0;70];
dir = 0; % degree
PlotCar(center, dir)

function PlotCar(center, dir)
    car_w = 20; % cm
    car_l = 15; % cm
    car_dir = [0;1];

    corner = [center(1)-car_w/2, center(2)-car_l/2; ...
              center(1)+car_w/2, center(2)-car_l/2; ...
              center(1)+car_w/2, center(2)+car_l/2; ...
              center(1)-car_w/2, center(2)+car_l/2];

    sensor = [center(1)-car_w/2+1, center(2); ...
              center(1)+car_w/2-1, center(2); ...
              center(1),   center(2)+car_l/2-1];

    R = [cosd(dir) -sind(dir); sind(dir) cosd(dir)]; % rotation matrix
    center = R*center;
    corner = (R*corner')';
    car_dir = R*car_dir;
    sensor = (R*sensor')';
    factor = 4;
    plot(center(1),center(2),'r+'); hold on
%     plot(sensor(:,1),sensor(:,2),'bo','MarkerSize',15,'MarkerFace','b');
%     plot(corner(:,1), corner(:,2),'bo','MarkerSize',10); hold on
    plot([corner(1,1) corner(2,1)],[corner(1,2) corner(2,2)],'LineWidth',1.5,'Color','k'); hold on
    plot([corner(2,1) corner(3,1)],[corner(2,2) corner(3,2)],'LineWidth',1.5,'Color','k'); hold on
    plot([corner(3,1) corner(4,1)],[corner(3,2) corner(4,2)],'LineWidth',1.5,'Color','k'); hold on
    plot([corner(1,1) corner(4,1)],[corner(1,2) corner(4,2)],'LineWidth',1.5,'Color','k'); hold on

    quiver(center(1),center(2),car_dir(1)*factor,car_dir(2)*factor,0,'LineWidth',2,'Color','g','MaxHeadSize',3);
    axis equal
end


