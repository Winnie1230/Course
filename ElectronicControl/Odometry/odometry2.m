clear all; close all; clc;

dt = 0.01;
% trapezoid parameters
seg1_t = 1;  seg2_t = 3;  seg3_t = 1;
seg_total = seg1_t + seg2_t + seg3_t;

%% command generator
% ----- command -----
vt = 2; % m/s
wt = 0.5; % rad/s
seg1_start = 0;
seg2_start = seg1_start + seg_total;
seg3_start = seg2_start + seg_total;
[t1,seg1_v] = GenerateTrapzoidCmd(seg1_start, seg1_t, seg2_t, seg3_t, vt, dt); % 0~7 sec
[t2,seg2_w] = GenerateTrapzoidCmd(seg2_start, seg1_t, seg2_t, seg3_t, wt, dt); % 7~14 sec
[t3,seg3_v] = GenerateTrapzoidCmd(seg3_start, seg1_t, seg2_t, seg3_t, vt, dt); % 0~7 sec
vcmd = [seg1_v; zeros(length(seg2_w),1) ; seg3_v];
wcmd = [zeros(length(seg1_v),1) ; seg2_w ; zeros(length(seg3_v),1)];
t = [t1 ; t2; t3];

%% coordinate calculate
theta = zeros(length(t),1); % car theta
x = zeros(length(t),1); % car position x
y = zeros(length(t),1); % car position y
for i = 2:length(t)
    theta(i) = theta(i-1) + wcmd(i)*dt;
    x(i) = x(i-1) + vcmd(i-1)*sin(theta(i-1))*dt;
    y(i) = y(i-1) + vcmd(i-1)*cos(theta(i-1))*dt;
end

%% plot animation
subplot(2,2,2);
h1 = animatedline('Color','b','LineWidth',1);
axis([0 seg_total*3 0 max(vcmd)]);
title('vcmd-t');
xlabel('t');
ylabel('vcmd');

subplot(2,2,4);
h2 = animatedline('Color','b','LineWidth',1);
axis([0 seg_total*3 0 max(wcmd)]);
title('wcmd-t plot');
xlabel('t(s)');
ylabel('wcmd');

subplot(2,2,3);
h3 = animatedline('Color','b','LineWidth',1);
axis([min(x) max(x) min(y) max(y)]);
title('x-y');
xlabel('x');
ylabel('y');

for i=1:length(t)
    addpoints(h1,t(i,1),vcmd(i,1));
    addpoints(h2,t(i,1),wcmd(i,1));
    addpoints(h3,x(i,1),y(i,1));
    drawnow;
end

%% function
function [t,cmd] = GenerateTrapzoidCmd(start_t, seg1_t, seg2_t, seg3_t, max_cmd, dt)
    % max_cmd: maximum cmd
    % segx_t: segment1 time duration(x=1,2,3)
    seg1_cmd = @(t)(max_cmd/seg1_t)*t;
    seg3_cmd = @(t) max_cmd - (max_cmd/seg3_t)*t;
    t1 = (0:dt:seg1_t)';
    t2 = (dt:dt:seg2_t)';
    t3 = (dt:dt:seg3_t)';
    cmd1 = seg1_cmd(t1);
    cmd2 = max_cmd*ones(length(t2),1);
    cmd3 = seg3_cmd(t3);
    t = start_t + [t1 ; seg1_t + t2 ; (seg1_t+seg2_t) + t3];
    cmd = [cmd1 ; cmd2 ; cmd3];
end