clear all; close all; clc;
% main function of SMD system
% compare different conditions(different w)
m = 1; c=5; k = 6;
t = 10; %sec
wn = (k/m)^(1/2); % natural frequency

% ----- PWM function parameters -----
A2 = [5 10 15 20]; 
alpha2 = 0.5;
w2 = 10;
% -----------------------------------

% ----- Experimental Group ------
for i = 1:length(A2)
    A = A2(i);
    
    % ----- state space representation -----
    f2 = @(t,X)[X(2);1/m*Force2(t,w2,alpha2,A)-k/m*X(1)-c/m*X(2)]; % PWM input
    % --------------------------------------
    [ts,xs] = ode45(f2,[0,t],[0;0]);
    
    subplot(2,1,1);
    plot(ts(:,1),xs(:,1),'DisplayName',['A=',num2str(A)]); hold on;
    title('x-t plot');
    xlabel('t');
    ylabel('x');
    legend;
    
    subplot(2,1,2);
    plot(ts(:,1),xs(:,2),'DisplayName',['A=',num2str(A)]); hold on;
    title('v-t plot')
    xlabel('t');
    ylabel('v');
    legend;
end


%% PWM function
function f = Force2(t,w,alpha,A2)
    tf = mod(t,2*pi/w);
    if tf <= alpha*2*pi/w
        f = A2;
    else
        f = 0;
    end
end