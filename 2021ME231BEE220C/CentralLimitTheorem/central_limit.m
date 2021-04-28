clear all; close all; clc;

iter = 1000;
r = unifrnd(-1,1,[iter,1]); % uniform distribution

k = 1024;
for i = 1 : iter
    rk(i,1) = 1/k^(1/2)*sum( r(randi([1,length(r)],k,1)));
end
nbins = 25;
[count,edges] = histcounts(rk, nbins);
dx = (edges(end) - edges(1))/length(count);
normalize_count = count/iter/dx;
figure
histogram('BinCounts', normalize_count, 'BinEdges', edges);
title(['Histogram of uniform distribution with sampling k = ',num2str(k)]); hold on;
xlabel(['x', num2str(k),'(sum over', num2str(k), ')']);
ylabel('Normalized count');

x = [-3:.1:3];
std(rk)
y = normpdf(x, 0, (1/3)^(1/2));
plot(x,y,'r'); hold off;
