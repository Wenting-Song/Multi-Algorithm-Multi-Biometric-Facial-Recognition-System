function [ CMC ] = calCMC( simMatrix, genuine, di, ti )
%CMC Summary of this function goes here
%   Detailed explanation goes here
CMC = zeros(100,1);

L = size(simMatrix,1);

if(di==0)
    Msorted = sort(simMatrix,2,'descend');
else
    Msorted = sort(simMatrix,2,'ascend');
end

for t = 1:100
    
    Mt = zeros(L,1);
    
    for i = 1:t
        Mt = Mt | (Msorted(:,i)==genuine);
    end
    
    Mtemp = Mt>0;
    CMC(t) = sum(Mtemp)/(L);
end

CMC=CMC.*100;

figure;
plot(CMC, 'g','linewidth',2.5);
title(ti);
ylabel('Rank-t Identification Rate (%)');
xlabel('Rank (t)');



end
