close all; clear all; clc;
% logistic map

%% variables 
Acoeff = [2]; % coefficient of logistic function
iter = 60; % number of iteration
init_x_min = 0.0001; % minimum value of initial x
init_x_max = 0.9999; % maximum value of initial x
initX = linspace(init_x_min,init_x_max,10);
y = @(A,x) (A*x*(1-x));
% ----- parameters of function FindFixPt -----
x0_arr = [0.1]; % initial value of x
iter1 = 50;
% --------------------------------------------

%% logistic map plot
for k = 1:size(Acoeff,2)
    A = Acoeff(k);
    for j = 1:size(initX,2)
        x(1,1) = initX(j);
        for i = 1:iter
            x(i+1,1) = y(A,x(i,1));
        end
        plot(x(1:iter,1),x(2:iter+1,1),'b*');
        
        % equilibrium point
        plot(0,0,'ro','markersize',7);
        plot(1-1/A,1-1/A,'ro','markersize',7);
        axis([0 1 0 1]);
        % axis equal;
        xlabel('x');
        ylabel('f(x)');
        hold on;
    end
    FindFixPt(A,x0_arr,iter1,y);
end

%% plot y=x
% x_n+1 = A*x_n*(1-x_n)
% y = x
xx = linspace(0,1,10);
yy = xx;
plot(xx,yy,'g','LineWidth',2);
hold off;

%% find fixed point
% find fixed point according to initial value of x
% One can use FindFixPt to check the fixed point is attractor or repeller
function FindFixPt(A,x0_arr,iter,y)
    % x0: initial value of x
    % iter: iteration of vertical lines and horizontal lines
    for j = 1:size(x0_arr,2)
        y0 = 0;
        x0 = x0_arr(j);
        for i = 1:iter
            y1 = y(A,x0);
            VertLine(x0,y0,y1);
            x1 = y1;
            HorizLine(y1,x0,x1);
            x0 = x1; y0 = y1;
        end
    end
end

% plot verticle line
function VertLine(x,yi,yf)
    %xi: start point  %xf: end point
    vy = linspace(yi,yf,10);
    vx = x*ones(1,size(vy,2));
    plot(vx,vy,'r','LineWidth',1);
end

% plot horizontal line
function HorizLine(y,xi,xf)
    %xi: start point  %xf: end point
    hx = linspace(xi,xf,10);
    hy = y*ones(1,size(hx,2));
    plot(hx,hy,'r','LineWidth',1);
end