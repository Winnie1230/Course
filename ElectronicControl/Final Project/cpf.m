clear all; close all; clc;

%% obstacle informations
obs_w = 25; % obstacle width
obs_l = 15; % obstacle length

% obstacle location
obs_loc = [5 20 ; -30 70; 5 120 ; -30 170];
obs_len = [obs_w obs_l; obs_w obs_l; obs_w obs_l; obs_w obs_l];
obs_center(:,1) = obs_loc(:,1)+obs_len(:,1)/2; 
obs_center(:,2) = obs_loc(:,2)+obs_len(:,2)/2;

w_coor = linspace(0,obs_w,obs_w+1);
l_coor = linspace(0,obs_l,obs_l+1);

% down/right/up/left divide equally and compute its coordinate
down = zeros(length(w_coor),2*length(obs_loc));
right = zeros(length(l_coor),2*length(obs_loc));
up = zeros(length(w_coor),2*length(obs_loc));
left = zeros(length(l_coor),2*length(obs_loc));

for i = 1:length(obs_loc)
    down(:,2*i-1:2*i) = obs_loc(i,:)+[w_coor',zeros(length(w_coor),1)];
    right(:,2*i-1:2*i) = obs_loc(i,:) + [obs_w*ones(length(l_coor),1),l_coor'];
    up(:,2*i-1:2*i) = obs_loc(i,:) + [w_coor',obs_l*ones(length(w_coor),1)];
    left(:,2*i-1:2*i) = obs_loc(i,:) + [zeros(length(l_coor),1),l_coor'];
end

%% Conventional Potential Field Method
P = [0;0]; % initial position
v = [0;0]; % initial velocity
m = 0.001; % mass
dt = 0.001; % sampling time
iter = 150;

target = [0;250];

katt = 5; % coefficient for attractive field
krep = 3; % coefficient for repulsive field
dmax = 30; % max tolarance distance between car and obstacle

h = figure('units','normalized','outerposition',[0 0 1 1]);
% filename = 'CPF.gif';
% del = 0.01; % time between animation frames

for j = 1:iter
    fatt = katt*(target-P)/norm(target-P);

    for i = 1:length(obs_loc)
        dir = (P-obs_center(i,:)');
    %     norm(dir)
        if(norm(dir) < dmax*2)
            if(dir(1)>=0 && dir(2)>=0) % first quadrant
                % up/right
                upframe = up(:,2*i-1:2*i);
                frep = krep*CalRep(upframe,dmax,P,krep);
                rightframe = right(:,2*i-1:2*i);
                frep = frep + CalRep(rightframe,dmax,P,krep);
            end
            if(dir(1)<=0 && dir(2)>=0) % second quadrant
                % up/left
                upframe = up(:,2*i-1:2*i);
                frep = krep*CalRep(upframe,dmax,P,krep);
                leftframe = left(:,2*i-1:2*i);
                frep = frep + CalRep(leftframe,dmax,P,krep);
            end
            if(dir(1)<=0 && dir(2)<=0) % third quadrant
                % down/left
                downframe = down(:,2*i-1:2*i);
                frep = krep*CalRep(downframe,dmax,P,krep);
                leftframe = left(:,2*i-1:2*i);
                frep = frep + CalRep(leftframe,dmax,P,krep);

            end
            if(dir(1)>=0 && dir(2)<=0) %forth quadrant
                % down/right
                downframe = down(:,2*i-1:2*i);
                frep = krep*CalRep(downframe,dmax,P,krep);
                rightframe = right(:,2*i-1:2*i);
                frep = frep + CalRep(rightframe,dmax,P,krep);
            end
        end
    end

    plot(P(1),P(2),'bo','MarkerSize',10,'LineWidth',4); hold on
    PlotQuiver(P,fatt,'g'); hold on
    PlotQuiver(P,frep,'r'); hold on
    PlotMap(obs_loc,obs_len); hold off
%     PlotQuiver(P,frep'+fatt,'b')
    axis equal
    axis([-50 50 -6 250])
    drawnow;
    
    F_total = fatt' + frep;
    a = F_total/m;
    % next timestep velocity and position
    v = v + a'*dt;
    P = P + v*dt + 1/2*a'*dt;
    
    if(norm(P-target)<1)
        break;
    end
    
    % ----- save to gif -----
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if j == 1
%         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
%     end
end

%% function
% Calculate down/right/up/left repulsive force
function frep = CalRep(frame,dmax,P,krep)
    dis = vecnorm(P - frame',2,1)';
    dir = (P - frame(dis < dmax,:)')';
    dir = dir./vecnorm(dir,2,2);
    dis = dis(dis < dmax);
    frep = krep*(sum(((1./dis) - 1/dmax).*dir,1));
end

function PlotQuiver(p1,dp,c)
%     dp = p2 - p1
    quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'Color',c,'MaxHeadSize',3);
end
