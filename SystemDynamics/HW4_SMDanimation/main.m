clear all; close all; clc;
%% animation
% ----- function parameters -----
A = 10; 
alpha = 0.6;
w = 3;
% -----------------------------------
ic = [4;2];
t = 5;
m = 1; c=5; k = 6;
SMD_Sineinput(m,c,k,A,w,t,ic);
SMD_PWMinput(m,c,k,A,alpha,w,t,ic);

m = 1; c=2; k = 2;
SMD_Sineinput(m,c,k,A,w,t,ic);
SMD_PWMinput(m,c,k,A,alpha,w,t,ic);

%% phase portrait
m = 1; c=5; k = 6;
p0 = -3:0.5:3; %initial condition range
PhasePortrait(m,c,k,p0);

m = 1; c = 2; k = 2;
PhasePortrait(m,c,k,p0);