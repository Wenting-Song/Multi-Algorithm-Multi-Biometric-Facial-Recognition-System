function [] = plotdist( Mimp, Mgen, lval, rval, t )
%PLOTDIST Summary of this function goes here
%   Detailed explanation goes here
Mimpdist = fitdist (Mimp, 'Normal');
Mgendist = fitdist (Mgen, 'Normal');

%plot(M1impdist);

x = lval:0.005:rval;

yimp = pdf(Mimpdist, x);
ygen = pdf(Mgendist, x);

figure();

h = plot(x, yimp, x, ygen);
title(t)
ylabel('Probability p(s)');
xlabel('Match Score (s)');
legend('Imposter Distribution', 'Genuine Distribution');
set(h(1),'linewidth',2.5);
set(h(2),'linewidth',2.5);

d = ((2^0.5)*abs(Mimpdist.mu - Mgendist.mu))/((Mimpdist.sigma^2 + Mgendist.sigma^2)^0.5)

end

