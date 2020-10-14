clear all; close all; clc;
% main function of SMD system
% use ode45 to solve solution
% compare different conditions(different w and different A)
m = 1; c=2; k = 2;
t = 4; %sec

% ----- sin function parameters -----
A1 = [5 10 15];
w1 = 10;
% -----------------------------------

for i = 1:length(A1)
    w = w1;
    A = A1(i);
    
    % ----- state space representation -----
    f1 = @(t,X)[X(2);1/m*Force1(t,w,A)-k/m*X(1)-c/m*X(2)]; % sin input
    % --------------------------------------
    [ts,xs] = ode45(f1,[0,t],[0;0]);
    
    % ----- calculate force -----
    for j = 1:size(ts,1)
        fs(j,1) = Force1(ts(j,1),w,A);
    end
    % ---------------------------
    
    subplot(3,1,1);
    plot(ts(:,1),fs(:,1),'DisplayName',['A=',num2str(A)]); hold on;
    title('f-t plot');
    xlabel('t');
    ylabel('f');
    legend;
    
    subplot(3,1,2);
    plot(ts(:,1),xs(:,1),'DisplayName',['A=',num2str(A)]); hold on;
    title('x-t plot');
    xlabel('t');
    ylabel('x');
    legend;
    
    subplot(3,1,3);
    plot(ts(:,1),xs(:,2),'DisplayName',['A=',num2str(A)]); hold on;
    title('v-t plot')
    xlabel('t');
    ylabel('v');
    legend;
    clear fs;
end

%% one pulse sine function
function f = Force1(t,w,A)
    if t <= 2*pi/w;
        f = A*sin(w*t);
    else
        f = 0;
    end
end