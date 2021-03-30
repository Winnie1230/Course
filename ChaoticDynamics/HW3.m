clear all; close all; clc;

%% variable
% a = 1.4; b = 0.3;
a = 0.5501; b = 0.88; c = -0.6001; d = 2.32001; e = 2.32; f = 0.043;
iter = 20000; % number of iteration
trunc = 500;

%% quadratic map
% x = @(a,xy)(1+xy(2)-a*xy(1)^2);
% y = @(b,xy)(b*xy(1));

x = @(a,b,c,xy)(xy(1)^2-a*xy(2)^2+b*xy(1)+c*xy(2));
y = @(d,e,f,xy)(d*xy(1)*xy(2)+e*xy(1)+f*xy(2));

init = [-0.5,-0.5];
xy_result(1,:) = init;

ap = 0.7 : 0.005 : 0.735;
% bp = 0.69 : 0.005 : 0.725;
% cp = -0.58 : 0.01 : -0.51;
% dp = 7.20 : 0.02 : 7.34;
% ep = 2.4 : 0.004 : 2.428;
% fp = 0.05 : 0.001 : 0.057;

for j = 1 : length(ap)
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    a = ap(j);
    for i = 2:iter
%         xy_result(i,:) = [x(a,xy_result(i-1,:)),y(b,xy_result(i-1,:))];
        xy_result(i,:) = [x(a,b,c,xy_result(i-1,:)),y(d,e,f,xy_result(i-1,:))];
    end
    plot(xy_result(trunc:iter,1),xy_result(trunc:iter,2),'b*');
    axis equal;
    title({['Quadratic Map'];['a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d),', e=',num2str(e),', f=',num2str(f)]});
    xlabel('x');
    ylabel('y');
end