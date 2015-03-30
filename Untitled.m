clear all;
close all;
clc;

iPredTauCoeffs=[1];
iPredRhoCoeffs=[1 -1];
iCorrTauCoeffs=[1 0];
iCorrRhoCoeffs=[1 -1];

wBetaK = iCorrRhoCoeffs(1);

wStabilityPolynom = [-1*poly2sym(iPredTauCoeffs*wBetaK,'z') poly2sym(wBetaK*iPredRhoCoeffs,'z')-poly2sym(iCorrTauCoeffs,'z') poly2sym(iCorrRhoCoeffs,'z')];

wStabilityRoots = roots(wStabilityPolynom);

wt=0:.01:2*pi;
x=exp(i*wt);

wRootsValues = cell(1,length(wStabilityRoots));
disp('Now computing value, long process, please stay tuned')
tic
for i=1:length(wStabilityRoots)
    wRootsValues{i} = subs(wStabilityRoots(i),x);
end
toc
disp('Values computed, now ploting.')



