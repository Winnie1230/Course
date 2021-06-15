clear all; close all; clc;

dt = 0.01;
% trapezoid parameters
v_seg1_t = 1;  v_seg2_t = 7;  v_seg3_t = 1;
w_seg1_t = 1;  w_seg2_t = 3;  w_seg3_t = 1;

%% command generator
% ----- command -----
vt = 1; % m/s
wt = 0.3; % rad/s
v_start = 0;
w_start = 2;
[tv,vcmd] = GenerateTrapzoidCmd(v_start, v_seg1_t, v_seg2_t, v_seg3_t, vt, dt);
[tw,seg_w] = GenerateTrapzoidCmd(w_start, w_seg1_t, w_seg2_t, w_seg3_t, wt, dt);
wcmd = zeros(length(tv),1);
wcmd((w_start/dt)+1: (w_start/dt)+1+length(seg_w)-1) = seg_w;
t = tv;
% subplot(2,1,1);
% plot(t,vcmd);
% subplot(2,1,2);
% plot(t,wcmd);

%% compute left, right motor command
d = 0.15; % distance between left and right wheel
r = 0.05; % radius of wheel
vcmd_r = vcmd + wcmd*d;
vcmd_l = vcmd - wcmd*d;

%% compute the car state using the measurement of left and right motor velocity
% now assume that the left and right motor can get to its ideal velocity at
% every time step
v_r = vcmd_r;
v_l = vcmd_l;

% we can compute the velocity and w of the car center(from measurement)
% but now use the ideal condition
% v_c = (v_r + v_l)/2;
% w_c = (v_r - v_l)/(2*d);

%% coordinate calculate
car_ideal_p = zeros(length(t),2); % car ideal position
theta_ideal = zeros(length(t),1); % car ideal theta
car_p_m = zeros(length(t),2);   % real car position
theta_m = zeros(length(t),1); % theta from measurement
wheel_r_m = zeros(length(t),2); % position of right wheel
wheel_l_m = zeros(length(t),2); % position of left wheel

wheel_r_m(1,:) = [d,0];
wheel_l_m(1,:) = [-d,0];

dis = zeros(length(t),1);

for i = 2:length(t)
    % ideal position
    theta_ideal(i) = theta_ideal(i-1) + wcmd(i-1)*dt;
    car_ideal_p(i,:) = car_ideal_p(i-1,:) + vcmd(i-1)*[-sin(theta_ideal(i-1)), cos(theta_ideal(i-1))]*dt;
    
    v_c = (v_r(i) + v_l(i))/2;
    w_c = (v_r(i) - v_l(i))/(2*d);
    
    theta_m(i) = theta_m(i-1) + w_c*dt;
    car_p_m(i,:) = car_p_m(i-1,:) + v_c*[-sin(theta_m(i-1)), cos(theta_m(i-1))]*dt;
    wheel_r_m(i,:) = car_p_m(i,:) + d*[cos(theta_m(i)), sin(theta_m(i))];
    wheel_l_m(i,:) = car_p_m(i,:) - d*[cos(theta_m(i)), sin(theta_m(i))];
    
    dis(i) = norm(wheel_r_m(i,:)-wheel_l_m(i,:));
end

%% static plot
% scatter(wheel_r_m(1:10:end,1),wheel_r_m(1:10:end,2)); hold on
% scatter(wheel_l_m(1:10:end,1),wheel_l_m(1:10:end,2)); hold on
% scatter(car_p_m(1:10:end,1),car_p_m(1:10:end,2)); hold on
% axis equal
figure('Name','Odometry','units','normalized','outerposition',[0 0 1 1])
subplot(2,2,2);
plot(t,vcmd,'Color','b','LineWidth',1);
title('vcmd-t');
xlabel('t[s]');
ylabel('vcmd[cm/s]');
% axis([0 t(end,1) 0 max(vcmd)+0.5]);

subplot(2,2,4);
plot(t,wcmd,'Color','b','LineWidth',1);
% axis([0 t(end,1) 0 max(wcmd)+0.1]);
title('wcmd-t plot');
xlabel('t[s]');
ylabel('wcmd[rad/s]');

subplot(2,2,[1,3]);
plot(car_p_m(:,1),car_p_m(:,2),'Color','b','LineWidth',1,'DisplayName','car center position'); hold on
plot(wheel_r_m(:,1),wheel_r_m(:,2),'Color','r','LineWidth',1,'DisplayName','right wheel position'); hold on
plot(wheel_l_m(:,1),wheel_l_m(:,2),'Color','g','LineWidth',1,'DisplayName','left wheel position'); hold off
axis([min(car_p_m(:,1))-1.5*d max(car_p_m(:,1))+1.5*d min(car_p_m(:,2)) max(car_p_m(:,2))]);
title('x-y');
xlabel('x');
ylabel('y');
legend
axis equal


%% plot animation
% figure('Animation');
% subplot(2,2,2);
% h1 = animatedline('Color','b','LineWidth',1);
% axis([0 t(end,1) 0 max(vcmd)]);
% title('vcmd-t');
% xlabel('t');
% ylabel('vcmd');
% 
% subplot(2,2,4);
% h2 = animatedline('Color','b','LineWidth',1);
% axis([0 t(end,1) 0 max(wcmd)]);
% title('wcmd-t plot');
% xlabel('t');
% ylabel('wcmd');
% 
% subplot(2,2,[1,3]);
% h3 = animatedline('Color','b','LineWidth',1);
% axis([min(car_p_m(:,1))-1.5*d max(car_p_m(:,1))+1.5*d min(car_p_m(:,2)) max(car_p_m(:,2))]);
% title('x-y');
% xlabel('x');
% ylabel('y');
% 
% subplot(2,2,[1,3]);
% h4 = animatedline('Color','r','LineWidth',1);
% axis([min(car_p_m(:,1))-1.5*d max(car_p_m(:,1))+1.5*d min(car_p_m(:,2)) max(car_p_m(:,2))]);
% title('x-y');
% xlabel('x');
% ylabel('y');
% 
% subplot(2,2,[1,3]);
% h5 = animatedline('Color','g','LineWidth',1);
% axis([min(car_p_m(:,1))-1.5*d max(car_p_m(:,1))+1.5*d min(car_p_m(:,2)) max(car_p_m(:,2))]);
% title('x-y');
% xlabel('x');
% ylabel('y');
% axis equal
% 
% for i=1:length(t)
%     addpoints(h1,t(i,1),vcmd(i,1));
%     addpoints(h2,t(i,1),wcmd(i,1));
%     addpoints(h3,car_p_m(i,1),car_p_m(i,2));
%     addpoints(h4,wheel_r_m(i,1),wheel_r_m(i,2));
%     addpoints(h5,wheel_l_m(i,1),wheel_l_m(i,2));
%     drawnow;
% end

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