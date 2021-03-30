clear all; close all; clc;
r = 1;
theta = (0:0.01:2*pi)';
lamda1 = 1;     ev1 = [1;1];
lamda2 = 0.01 ; ev2 = [-1;1];

x = r*cos(theta); y = r*sin(theta);
coor0 = [x,y];
M = [ev1,ev2]*[lamda1 0 ; 0 lamda2]*inv([ev1,ev2]);

coor1 = coor0*M;

l1 = (1+(ev1(1)/ev1(2))^2)^(1/2);
% ang2 = acos(dot(ev2,[1;0])/norm(ev2));
% l2 = 0.2;
% pf2 = [r*cos(ang2),r*sin(ang2)];
% pi2 = pf2 + l2*ev2/norm(ev2);
% dp2 = pf2 - pi2;
figure('units','normalized','outerposition',[0 0 1 1])
plot(coor1(:,1),coor1(:,2),'Color',[0 0.4470 0.7410],'LineWidth',2); hold on;
% plot(x,y,'Color',[0 0.4470 0.7410],'LineWidth',2); hold on;
plot(0,0,'Marker','+','Color','black','MarkerSize',20,'LineWidth',2); hold on;
plot([-l1*ev1(1)/norm(ev1),l1*ev1(1)/norm(ev1)],[-l1*ev1(2)/norm(ev1),l1*ev1(2)/norm(ev1)],'r','LineWidth',2,'LineStyle','--'); hold on;
% plot([-l2*ev2(1)/norm(ev2),l2*ev2(1)/norm(ev2)],[-l2*ev2(2)/norm(ev2),l2*ev2(2)/norm(ev2)],'Color',[0.8500 0.3250 0.0980],'LineWidth',2,'LineStyle','--');
% quiver(pi2,pf2,dp2(1),dp2(2),0);
axis([-1 1 -1 1]);
xlabel('x');
ylabel('y');
axis equal tight;



