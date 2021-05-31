clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework1
% Generate desired PDF of continuous random variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda = 1000;
f = @(t) 2*t/(lamda^2).*exp(-t.^2/(lamda^2)); % pdf

t_range = 3000;
% vpasolve(integral(f,0,T)==0.5,T)
% q = integral(f,0,422.9799)
t1 = linspace(0,t_range,t_range+1)';
cdf = zeros(length(t1),1);
pdf = f(t1);
for i = 1:length(t1)
    cdf(i,1) = integral(f,0,t1(i));
end
% subplot(1,2,1); plot(t1, pdf);
% title('PDF');
% xlabel('PDF f(t)');
% ylabel('Time to failure t[days]')
% subplot(1,2,2); plot(t1, cdf);
% title('CDF');
% xlabel('t[days]');
% ylabel('probability');
sample_num = 1000000; % sample numbers

r = unifrnd(0,1,[sample_num,1]); % uniform distribution
% ----- probability of failure r before T -----
T = lamda * (-log(1-r)).^(1/2);
% ---------------------------------------------

% ----- information of histogram -----
nbins = t_range / 10; % number of bins in histogram
[count,edges] = histcounts(T, nbins);
dx = (edges(end) - edges(1))/length(count);
normalize_count = count/sample_num/dx; % normalize count
% ------------------------------------

histogram('BinCounts', normalize_count, 'BinEdges', edges, 'FaceColor', 'b'); hold on;
plot(t1, pdf,'r','LineWidth',1.5);
title('PDF');
xlabel('PDF f(t)');
ylabel('Time to failure t[days]'); hold off

tm = 1;


