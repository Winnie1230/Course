close all; clear all; clc;
%% parameters
% unit: length(m), mass(kg)
H = 1.9; % height of railing
l1 = 1.60; % height
l2 = 0.3; % length of down part of arm
l3 = 0.3; % length of up part of arm
L = 0.3;
m = 50; % weight
g = 9.8; %m/s^2

d0 = (L^2 + (H-l1)^2)^(1/2); % origin length of spring
theta = 0; % rad

[alpha,beta,phi,d] = CalLengthAngle(l1,l2,l3,H,L,theta);

%% plot
plot([0,-l1*sin(theta)],[0,l1*cos(theta)],'k','LineWidth',2); hold on;
plot([-l1*sin(theta),l2*cos(alpha)-l1*sin(theta)],[l1*cos(theta),l1*cos(theta)+l2*(alpha)],'k','LineWidth',2); hold on;
plot([l2*cos(alpha)-l1*sin(theta),L],[l1*cos(theta)+l2*(alpha),H],'k','LineWidth',2); hold on;

plot(0,0,'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;
plot(-l1*sin(theta),l1*cos(theta),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;
plot(l2*cos(alpha)-l1*sin(theta),l1*cos(theta)+l2*(alpha),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;
plot(L,H,'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;

axis equal;
xlim([-0.3 L+0.2]);
ylim([0 H]);

%% function
function [alpha,beta,phi,d] = CalLengthAngle(l1,l2,l3,H,L,theta)
    % ----- constraint(max theta) -----
    t1 = asin((H^2+L^2+l1^2-l2^2-l3^2-2*l2*l3)/(2*l1*(H^2+L^2)^(1/2)));
    t2 = asin(H/(H^2+L^2)^(1/2));
    theta_max = t2 - t1;
    if(theta > theta_max)
        theta = theta_max;
    end
    % ---------------------------------
    d = ((H-l1*cos(theta))^2+(L+l1*sin(theta))^2)^(1/2); % spring length  
    beta = atan2(H-l1*cos(theta),L+l1*sin(theta));
    theta1 = acos((d^2+l2^2-l3^2)/(2*d*l2));
    theta2 = acos((d^2+l3^2-l2^2)/(2*d*l3));
    alpha = beta - theta1;
    phi = (pi/2 - beta) - theta2;
end
