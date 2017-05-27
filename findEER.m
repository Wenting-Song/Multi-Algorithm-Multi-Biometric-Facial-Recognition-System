function [ EER, thr ] = findEER( FMR, FNMR )
%EER Summary of this function goes here
%   Detailed explanation goes here

% Find EER
d = (2.*((FMR-FNMR).^2)).^0.5;
[~,pos] = min(d);
EER = FMR(pos);
T=0:0.0001:1;
thr = T(pos);

end

