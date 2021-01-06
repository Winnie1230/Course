clear all; close all; clc;
%% phase portrait
% assume f(t) = 0, draw x-x_dot phase portrait using lots of initial conditions
m = 1; c = 2; k = 2;
% ----- state space representation -----
f2 = @(t,X)[X(2);-k/m*X(1)-c/m*X(2)]; % sin input
% --------------------------------------
figure('Name','SMD system_phase portrait','units','normalized','outerposition',[0 0 1 1])
p0 = -3:0.5:3; %initial condition range
for i = 1:length(p0)
    [ts,ys] = ode45(f2,[0 20],[p0(i);(-1)^i*p0(i)]);
    plot(ys(:,1),ys(:,2));  hold on
    plot(ys(1,1),ys(1,2),'bo'); hold on % starting point
    plot(ys(end,1),ys(end,2),'kx','MarkerSize',15,'LineWidth',4); hold on % ending point
    grid on
    axis equal tight;
end
hold off;