function [OUTPUTT]= ISS_NNN(SPOTIS)
%%%%ISS NNN creates the Nearest Neighbor Network

NUMBER=[];
OWN=[];
NEAR=[];
for elem=1:size(SPOTIS.pos,1)
    disp(elem)
    DISTANCES=sqrt((SPOTIS.pos(elem,1)-SPOTIS.pos(:,1)).^2 + (SPOTIS.pos(elem,2)-SPOTIS.pos(:,2)).^2);
    COMB=table(DISTANCES,SPOTIS.name);
    B = sortrows(COMB,1);
    for i =2:4
        NUMBER=[NUMBER,elem];
        OWN=[OWN,SPOTIS.name(elem)];
        NEAR=[NEAR,B{i,2}];
    end   
end

COEXPREL=table(NUMBER',OWN',NEAR');

[coexp,chi2,p,label]=crosstab(COEXPREL.Var2,COEXPREL.Var3);

coexp=array2table(coexp);
coexp.Properties.RowNames=label(:,1);
coexp.Properties.VariableNames=label(:,2);
coexp = sortrows( coexp,'RowNames');
coexp=coexp(:,sort(coexp.Properties.VariableNames));

%SIMULATIONS

NUMBERSIM=[];
OWNSIM=[];
NEARSIM=[];

SIMSPOTIS=SPOTIS;
SIMSPOTIS.name=SIMSPOTIS.name(randperm(length(SIMSPOTIS.pos)));
for elem=1:size(SIMSPOTIS.pos,1)
    disp(elem)
    DISTANCES=sqrt((SIMSPOTIS.pos(elem,1)-SIMSPOTIS.pos(:,1)).^2 + (SIMSPOTIS.pos(elem,2)-SIMSPOTIS.pos(:,2)).^2);
    COMB=table(DISTANCES,SIMSPOTIS.name);
    B = sortrows(COMB,1);
    for i =2:4
        NUMBERSIM=[NUMBERSIM,elem];
        OWNSIM=[OWNSIM,SIMSPOTIS.name(elem)];
        NEARSIM=[NEARSIM,B{i,2}];
    end   
end


COEXPRELSIM=table(NUMBERSIM',OWNSIM',NEARSIM');

[coexpSIM,chi2,p,label2]=crosstab(COEXPRELSIM.Var2,COEXPRELSIM.Var3);


coexpSIM=array2table(coexpSIM);

coexpSIM = sortrows( coexpSIM,'RowNames');
coexpSIM=coexpSIM(:,sort(coexpSIM.Properties.VariableNames));

ENRICH=(table2array(coexp)-table2array(coexpSIM))./((table2array(coexp)+table2array(coexpSIM))./2);


figure
heatmap(coexp.Properties.VariableNames,coexp.Properties.VariableNames,ENRICH);
set(gca,'Colormap',jet(30))

OUTPUTT=array2table(ENRICH);
OUTPUTT.Properties.RowNames=coexpSIM.Properties.VariableNames;
OUTPUTT.Properties.VariableNames=coexpSIM.Properties.VariableNames;






end