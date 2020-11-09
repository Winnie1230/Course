close all; clear all; clc;

% w = 10;
% t = 0:0.01:5;
% 
% f = zeros(1,length(t));
% for i = 1:2:5
%     f = f + 2/(i*pi)*sin(i*w.*t);
% end
% 
% plot(t,f);

H = tf([1],[1 5 6]);
bode(H)
fb = bandwidth(H)