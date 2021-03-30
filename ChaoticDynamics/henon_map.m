clear all; close all; clc;

%% variable
% a = 1.4; b = 0.3;

a = 0.02;
iter = 20000; % number of iteration
trunc = 500;

%% quadratic map
% x = @(a,xy)(1+xy(2)-a*xy(1)^2);
% y = @(b,xy)(b*xy(1));

x = @(xy)(-0.3-0.13*xy(1)+0.1*xy(2)-0.8*xy(1)*xy(2)-0.5*xy(1)^2+0.4*xy(2)^2);
y = @(a,xy)(-1+1.1*xy(1)-0.3*xy(2)-1.6*xy(1)*xy(2)+1.3*xy(1)^2+a*xy(2)^2);

init = [0,0];
xy_result(1,:) = init;

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 2:iter
    xy_result(i,:) = [x(xy_result(i-1,:)),y(a,xy_result(i-1,:))];
end
plot(xy_result(trunc:iter,1),xy_result(trunc:iter,2),'b*');
% axis equal tight;
% title({['Quadratic Map'];['a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)]});
xlabel('x');
ylabel('y');
grid on;
