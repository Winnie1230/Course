clear all; close all; clc;
% pdf of a uniform die
% x: discrete random variable(finite outcome space)
% y: probability of every outcome
x = [1,2,3,4,5,6];
y = 1/length(x)*ones(length(x),1);

cdf = zeros(1,length(x));
for i = 1:length(x)
    cdf(1,i) = sum(y(1:i));
end

figure();
subplot(1,2,1);
plot(x,y,'bo','MarkerSize',7); hold on
for i = 1:length(x)
    plot([x(i),x(i)],[0 y(i)],'LineWidth',1.5,'Color','r'); hold on;
end
axis([1 6 0 0.3]);
title('pdf of a uniform die')
xlabel('die count');
ylabel('probability'); hold off;

subplot(1,2,2);
plot(x,cdf,'bo');hold on;
plot(x(2:end),cdf(1:end-1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b'); hold on
for i = 1:length(x)-1
    plot([x(i),x(i+1)],[cdf(i),cdf(i)],'LineWidth',1.5,'Color','r'); hold on;
end
axis([1 6 0 1]);
title('cdf of a uniform die')
xlabel('die count');
ylabel('cdf');