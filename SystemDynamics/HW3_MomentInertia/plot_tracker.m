clear all; close all; clc;

xlsFile1 = 'tracker_g.xlsx';
xlsFile2 = 'tracker_b.xlsx';
data1 = xlsread(xlsFile1);
data2= xlsread(xlsFile2);
data2 = rmmissing(data2);

figure
[period1,period_avg1,period_std1] = PlotTracker(data1);
figure
[period2,period_avg2,period_std2] = PlotTracker(data2);

function [period,period_avg,period_std] = PlotTracker(data)
    t = seconds(linspace(0,max(data(:,1)),length(data(:,1))));
    min_pos = min(data(:,4));
    TF = islocalmin(data(:,4),'MinSeparation',seconds(0.6),'SamplePoints',t);

    plot(data(:,1),data(:,4)-min_pos,'b-','DisplayName','tracker data'); hold on;
    plot(data(TF,1),data(TF,4)-min_pos,'r*','DisplayName','local minimum');
    title('Tracker data plot');
    xlabel('t(s)');
    ylabel('theta(deg)');
    legend;
    grid on;
    axis tight;

    local_min = data(TF,1);
    period = zeros(length(local_min)-1,1);

    for i = 1 : length(period)
        period(i) = local_min(i+1)-local_min(i);
    end

    period_std = std(period);
    period_avg = mean(period);
end