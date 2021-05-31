clear all; close all; clc;
cr = 250;
cm = 50;
lamda = 1000;
c = cr^2*pi/(cm^2*lamda^2);
% sval = fzero(@(t)c*t^2+exp(-t^2/(lamda^2))-1, 1);
sval = fzero(@(t)2*cr*pi^(1/2)*t^3 - cm*lamda^3*(1-exp(-t^2/(lamda^2)))^(1/2)*(1+2*t^2/(lamda^2)), 500);
% v = 1; % Just an example value
% f = @(t) 2*t/(lamda^2).*exp(-t.^2/(lamda^2)); % pdf
% sval = fzero(@(s)integral(f,t,0,s)-v, 1)