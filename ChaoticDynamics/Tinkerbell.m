clear all; close all; clc;

%% variable
% a = 1.4; b = 0.3;
a = 0.9; b = -0.6013; c = 2.0; d = 0.5;
iter = 20000; % number of iteration
trunc = 500;

%% quadratic map
x = @(a,b,xy)(xy(1)^2-xy(2)^2+a*xy(1)+b*xy(2));
y = @(c,d,xy)(2*xy(1)*xy(2)+c*xy(1)+d*xy(2));

init = [-0.72,-0.64];
xy_result(1,:) = init;

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 2:iter
    xy_result(i,:) = [x(a,b,xy_result(i-1,:)),y(c,d,xy_result(i-1,:))];
end
plot(xy_result(trunc:iter,1),xy_result(trunc:iter,2),'b*');
axis equal tight;
% title({['Quadratic Map'];['a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)]});
xlabel('x');
ylabel('y');
grid on;
