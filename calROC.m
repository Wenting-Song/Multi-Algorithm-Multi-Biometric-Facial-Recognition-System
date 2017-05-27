function [ FMR, FNMR ] = calROC( imposter, genuine,thr, di, ti )
%CALROC Summary of this function goes here
%   Detailed explanation goes here

T = 0:0.0001:thr;
FNMR = [];

% Calculate FNMR

for i = T
    if(di == 0)
        FNMR = [FNMR ; (sum(genuine<i) / length(genuine))];
    else
        FNMR = [FNMR ; (sum(genuine>i) / length(genuine))];
    end
end

GMR = 1-FNMR;

% calculate FMR

FMR = [];

for i = T
    if(di == 0)
        FMR = [FMR ; (sum(imposter>i)/(length(imposter)))];
    else
        FMR = [FMR ; (sum(imposter<i)/(length(imposter)))];
    end
end

FMR = FMR.*100;
FNMR = FNMR.*100;

figure;
plot(FMR,FNMR, 1:100, 1:100,'--g','linewidth', 2.5);
legend('FAR vs FRR', 'y=x');
title(ti);
ylabel('False Reject(Non-Match) Rate (%)');
xlabel('False Accept(Match) Rate (%)');
end

