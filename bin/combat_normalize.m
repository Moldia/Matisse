% Normalize cells using Combat
% Sergio Marco, 2020
function[CELLS] =combat_normalize(BINS)
CELLIS=BINS;
%hist(max(BINS.expressionmatrix'),30)
CELLICS=CELLIS;
samplec=unique(CELLIS.sample);
for ss=1:size(samplec);
    namec=samplec(ss);
    CELLICS.sample(ismember(CELLIS.sample,namec))={ss};
end
SAMPLE=cell2mat(CELLICS.sample);
bayesdata=combat(CELLIS.expressionmatrix',SAMPLE,[] , 1);
CELLIS.expressionmatrix=bayesdata';
CELLS=CELLIS;
end