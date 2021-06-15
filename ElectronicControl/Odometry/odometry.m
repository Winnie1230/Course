clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 1
% divide velocity command into three segment with zero angular velocity
% segment 1: constant acceleration
% segment 2: constant velociyt
% segment 3: constant acceleratoin
% assume initial velocity = 0 and initial angle = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.01;
seg1_t = 1; % segment1 time duration
seg2_t = 7; % segment2 time duration
seg3_t = 1; % segment3 time duration

% ----- command -----
vt = 5; % rad/s

t1 = (0:dt:seg1_t)';
t2 = (0:dt:seg2_t)';
t3 = (0:dt:seg3_t)';

seg1_v = @(t)(vt/seg1_t)*t;
seg3_v = @(t) vt - (vt/seg3_t)*t;
v1 = seg1_v(t1);
v2 = vt*ones(length(t2),1);
v3 = seg3_v(t3);

t = [t1 ; seg1_t + t2 ; (seg1_t+seg2_t) + t3];
v = [v1 ; v2 ; v3];
w = zeros(length(t),1);

% plot(t,v);
% axis([0 (seg1_t+seg2_t+seg3_t) 0 vt+0.5]);
% title('velocity command');
% xlabel('t');
% ylabel('v');

theta = zeros(length(t),1); % car theta
x = zeros(length(t),1); % car position x
y = zeros(length(t),1); % car position y

for i = 2:length(t)
    theta(i) = theta(i-1) + w(i)*dt;
    x(i) = x(i-1) + v(i-1)*sin(theta(i-1))*dt;
    y(i) = y(i-1) + v(i-1)*cos(theta(i-1))*dt;
end

plot(x,y,'b-','LineWidth',1.5);

