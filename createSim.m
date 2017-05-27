function [ simMatrix, imposter, genuine ] = createSim( probe, gallery )
%CREATESIM Summary of this function goes here
%   Detailed explanation goes here

simMatrix = zeros(size(probe,2), size(gallery,2));
genuine = [];
imposter = [];

for i = 1:size(probe,2)
    for j = 1:size(gallery,2)
        simMatrix(i,j) = corr2(probe(:,i),gallery(:,j));
        if (ceil(i/2)==j)
            genuine = [genuine; simMatrix(i,j)];
        else
            imposter = [imposter; simMatrix(i,j)];
        end
    end
end
end

