clear all; close all; clc;
% Generate desired PDF of discrete random variable

lamda = 1000;
f = @(t) 2*t/(lamda^2).*exp(-t.^2/(lamda^2)); % pdf

t_range = 3000;
t1 = linspace(0,t_range,t_range+1)';
pdf = f(t1);
cdf = zeros(length(t1),1);

for i = 1:length(t1)
    cdf(i,1) = sum(pdf(1:i,1));
end

subplot(1,2,1);
plot(t1,pdf);
title('pdf');
xlabel('t');
ylabel('pdf f(t)');

subplot(1,2,2);
plot(t1,cdf);
title('cdf');
xlabel('t');
ylabel('probability');


sample_num = 1000000; % sample numbers
r = unifrnd(0,1,[sample_num,1]); % uniform distribution

% ----- probability of failure r before T -----
T = zeros(sample_num,1);
for i = 1:sample_num
    [row,col] = find(cdf < r(i));
    T(i,1) = t1(max(row),1);
end
% ---------------------------------------------

% ----- information of histogram -----
nbins = t_range / 10; % number of bins in histogram
[count,edges] = histcounts(T, nbins);
dx = (edges(end) - edges(1))/length(count);
normalize_count = count/sample_num/dx; % normalize count
% ------------------------------------

% histogram('BinCounts', normalize_count, 'BinEdges', edges, 'FaceColor', 'b'); hold on;
% plot(t1, pdf,'r','LineWidth',1.5);
% title('PDF');
% xlabel('PDF f(t)');
% ylabel('Time to failure t[days]'); hold off

tm = 1;   % preventative maintenance
cm = 50;  % preventative maintenance costs
cr = 250; % repairing cost

[row_below_index, col_below] = find(T <= tm & T > 0);
[row_beyond, col_beyond] = find(T > tm);

% ----- average running cost -----
Rr = cr./T(row_below_index,1);
Rm = (cm/tm)*ones(length(row_beyond),1);
R = mean([Rr;Rm]);

