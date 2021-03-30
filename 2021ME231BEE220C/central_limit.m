clear all; close all; clc;

iter = 500;
r = unifrnd(-1,1,[iter,1]);
% figure;
% h = histogram(r);
% title('Histogram of uniform distribution');

k = 64;
for i = 1 : iter
    rk(i,1) = 1/k^(1/2)*sum( r(randi([1,length(r)],k,1)));
end
figure
histogram(rk,'Normalization','probability');
title(['Histogram of uniform distribution with sampling k = ',num2str(k)]); hold on;
% h = histogram(r,'Normalization','probability'); hold off;

x = [-2:.1:2];
y = normpdf(x,mean(rk),std(rk));
plot(x,y,'r'); hold off;


% 
% gridsize = 0.5;
% unirange = 1;
% count = zeros(1,unirange/gridsize*2);
% iter = 10;
% 
% for i = 1 : iter
%     r = unifrnd(-1,1);
%     if mod(r,gridsize) ~= 0
%         count(1,fix(r/gridsize)+unirange/gridsize) = count(1,fix(r/gridsize)+unirange/gridsize) + 1;
%     else
%         count(1,fix(r/gridsize)+unirange/gridsize+1) = count(1,fix(r/gridsize)+unirange/gridsize+1) + 1;
%     end
%     
% end

