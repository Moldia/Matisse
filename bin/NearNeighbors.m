function [OUTPUTT]= Colocalization(CELLS,varargin)
%%%%ISS NNN creates the Nearest Neighbor Network

if nargin==1
imagesel==1;
else
 imagelse==varargin{2}
end

allsamples=unique(CELLS.sample);
SELECTED=ismember(CELLS.sample,allsamples(imagesel));

SELECTIONS=CELLS;
SELECTIONS.spotname=SELECTIONS.spotname(SELECTED,:);
SELECTIONS.location=SELECTIONS.location(SELECTED,:);
SELECTIONS.sample=SELECTIONS.sample(SELECTED,:);
if strcmp(class(SELECTIONS.images),'cell')
SELECTIONS.images=SELECTIONS.images{imagesel};
else
    SELECTIONS.images=SELECTIONS.images(imagesel);
end
SELECTIONS.scale=SELECTIONS.scale(imagesel);


SPOTIS.pos=SELECTIONS.location;
SPOTIS.name=SELECTIONS.spotname;


NUMBER=[];
OWN=[];
NEAR=[];
for elem=1:size(SPOTIS.pos,1)
    DISTANCES=sqrt((SPOTIS.pos(elem,1)-SPOTIS.pos(:,1)).^2 + (SPOTIS.pos(elem,2)-SPOTIS.pos(:,2)).^2);
    COMB=table(DISTANCES,SPOTIS.name);
    B = sortrows(COMB,1);
    for i =1:SELECTIONS.Colocalization_number_neighbors
        NUMBER=[NUMBER,elem];
        OWN=[OWN,SPOTIS.name(elem)];
        NEAR=[NEAR,B{i,2}];
    end   
end

COEXPREL=table(NUMBER',OWN',NEAR');

[coexp,chi2,p,label]=crosstab(COEXPREL.Var2,COEXPREL.Var3);

coexp=array2table(coexp);
coexp.Properties.RowNames=label(:,1);
coexp.Properties.VariableNames=label(1:size(coexp,2),2);
unadded=label(~ismember(label(1:size(coexp,1),1),label(1:size(coexp,2),2)),1);
for adding=1:size(unadded,1)
    coexp.Variable=zeros(size(coexp,1),1);
    coexp.Properties.VariableNames(end)=unadded(adding);
end

coexp = sortrows( coexp,'RowNames');
coexp=coexp(:,sort(coexp.Properties.VariableNames));

disp('Starting simulation');
%SIMULATION
coexpSIMTOT=zeros(size(coexp,1),size(coexp,2),100);
for total=1:100
disp(total);
NUMBERSIM=[];
OWNSIM=[];
NEARSIM=[];

SIMSPOTIS=SPOTIS;
SIMSPOTIS.name=SIMSPOTIS.name(randperm(length(SIMSPOTIS.pos)));
for elem=1:size(SIMSPOTIS.pos,1)
    DISTANCES=sqrt((SIMSPOTIS.pos(elem,1)-SIMSPOTIS.pos(:,1)).^2 + (SIMSPOTIS.pos(elem,2)-SIMSPOTIS.pos(:,2)).^2);
    COMB=table(DISTANCES,SIMSPOTIS.name);
    B = sortrows(COMB,1);
    for i =1:Colocalization_number_neighbors
        NUMBERSIM=[NUMBERSIM,elem];
        OWNSIM=[OWNSIM,SIMSPOTIS.name(elem)];
        NEARSIM=[NEARSIM,B{i,2}];
    end   
end


COEXPRELSIM=table(NUMBERSIM',OWNSIM',NEARSIM');

[coexpSIM,chi2,p,label2]=crosstab(COEXPRELSIM.Var2,COEXPRELSIM.Var3);


coexpSIM1=array2table(coexpSIM);
coexpSIM1.Properties.RowNames=label2(1:size(coexpSIM,1),1);
coexpSIM1.Properties.VariableNames=label2(1:size(coexpSIM,2),2);
coexpSIM1 = sortrows( coexpSIM1,'RowNames');
coexpSIM1=coexpSIM1(:,sort(coexpSIM1.Properties.VariableNames));
coexpSIMTOT(:,:,total)=table2array(coexpSIM1);
end

%%SIMULATIONS

% NUMBERSIM=[];
% OWNSIM=[];
% NEARSIM=[];
% 
% SIMSPOTIS=SPOTIS;
% SIMSPOTIS.name=SIMSPOTIS.name(randperm(length(SIMSPOTIS.pos)));
% for elem=1:size(SIMSPOTIS.pos,1)
%     disp(elem)
%     DISTANCES=sqrt((SIMSPOTIS.pos(elem,1)-SIMSPOTIS.pos(:,1)).^2 + (SIMSPOTIS.pos(elem,2)-SIMSPOTIS.pos(:,2)).^2);
%     COMB=table(DISTANCES,SIMSPOTIS.name);
%     B = sortrows(COMB,1);
%     for i =2:5
%         NUMBERSIM=[NUMBERSIM,elem];
%         OWNSIM=[OWNSIM,SIMSPOTIS.name(elem)];
%         NEARSIM=[NEARSIM,B{i,2}];
%     end   
% end
% 
% 
% COEXPRELSIM=table(NUMBERSIM',OWNSIM',NEARSIM');
% 
% [coexpSIM2,chi2,p,label2]=crosstab(COEXPRELSIM.Var2,COEXPRELSIM.Var3);
% 
% 
% coexpSIM2=array2table(coexpSIM2);
% coexpSIM2.Properties.RowNames=label2(1:size(coexpSIM2,1),1);
% coexpSIM2.Properties.VariableNames=label2(1:size(coexpSIM2,2),2);
% coexpSIM2 = sortrows( coexpSIM2,'RowNames');
% coexpSIM2=coexpSIM(:,sort(coexpSIM2.Properties.VariableNames));
% 
% %%SIMULATIONS
% 
% NUMBERSIM=[];
% OWNSIM=[];
% NEARSIM=[];
% 
% SIMSPOTIS=SPOTIS;
% SIMSPOTIS.name=SIMSPOTIS.name(randperm(length(SIMSPOTIS.pos)));
% for elem=1:size(SIMSPOTIS.pos,1)
%     disp(elem)
%     DISTANCES=sqrt((SIMSPOTIS.pos(elem,1)-SIMSPOTIS.pos(:,1)).^2 + (SIMSPOTIS.pos(elem,2)-SIMSPOTIS.pos(:,2)).^2);
%     COMB=table(DISTANCES,SIMSPOTIS.name);
%     B = sortrows(COMB,1);
%     for i =2:5
%         NUMBERSIM=[NUMBERSIM,elem];
%         OWNSIM=[OWNSIM,SIMSPOTIS.name(elem)];
%         NEARSIM=[NEARSIM,B{i,2}];
%     end   
% end
% 
% 
% COEXPRELSIM=table(NUMBERSIM',OWNSIM',NEARSIM');
% 
% [coexpSIM3,chi2,p,label2]=crosstab(COEXPRELSIM.Var2,COEXPRELSIM.Var3);
% 
% 
% coexpSIM3=array2table(coexpSIM3);
% coexpSIM3.Properties.RowNames=label2(1:size(coexpSIM3,1),1);
% coexpSIM3.Properties.VariableNames=label2(1:size(coexpSIM3,2),2);
% coexpSIM3 = sortrows( coexpSIM3,'RowNames');
% coexpSIM3=coexpSIM(:,sort(coexpSIM3.Properties.VariableNames));
% 
% 
% 
% %coexpSIM=(table2array(coexpSIM)+table2array(coexpSIM2)+table2array(coexpSIM3))./3;
% %coexpSIM=array2table(coexpSIM);
% %coexpSIM.Properties.RowNames=coexpSIM2.Properties.RowNames;
% %coexpSIM.Properties.VariableNames=coexpSIM2.Properties.VariableNames;
% 
% %CAR3=zeros(size(coexpSIM,1),1);
% %coexpSIM.Car3=CAR3;
% %coexpSIM=coexpSIM(:,sort(coexpSIM.Properties.VariableNames));
% A=[]
% A(:,:,1)=table2array(coexpSIM);
% A(:,:,2)=table2array(coexpSIM2);
% A(:,:,3)=table2array(coexpSIM3);
% A=coexpSIMTOT;
% pvals=zeros(size(coexp,1),size(coexp,2));
% for x=1:size(coexp,1)
%     for y=1:size(coexp,2)
%        same=A(x,y,:);
%       [sig,pval]=ttest(same,table2array(coexp(x,y)),'Tail','right');
%       pvals(x,y)=pval; 
%     end
% end
% 
% 
% %figure
% %clustergram(pvals,'Columnlabels',coexp.Properties.VariableNames,'Rowlabels',coexp.Properties.VariableNames)
% 
% %LEFT TEST
% 
% A=coexpSIMTOT;
% pvals2=zeros(size(coexp,1),size(coexp,2));
% for x=1:size(coexp,1)
%     for y=1:size(coexp,2)
%        same=A(x,y,:);
%       [sig,pval]=ttest(same,table2array(coexp(x,y)),'Tail','left');
%       pvals2(x,y)=pval; 
%     end
% end
% 
% figure
% subplot(1,2,1)
% he=heatmap(coexp.Properties.VariableNames,coexp.Properties.VariableNames,pvals);
% he.Colormap=hot;
% 
% subplot(1,2,2)
% h=heatmap(coexp.Properties.VariableNames,coexp.Properties.VariableNames,pvals2);
% h.Colormap=cool;
% 
% 
% SUMMARY=zeros(size(pvals,1),size(pvals,2));
% for xi=1:size(pvals,1)
%     for yi=1:size(pvals,2)
%     if (1-pvals(xi,yi))>(1-pvals2(xi,yi))
%         SUMMARY(xi,yi)=(1-pvals(xi,yi)).^100;
%     else
%         SUMMARY(xi,yi)=-((1-pvals2(xi,yi)).^100);
%         
%     end  
%     end
% end
% 
% figure
% h=heatmap(coexp.Properties.VariableNames,coexp.Properties.VariableNames,SUMMARY);
% h.Colormap=redbluecmap(30);
% 



A=coexpSIMTOT;
DIFF=zeros(size(coexp,1),size(coexp,2));
for x=1:size(coexp,1)
    for y=1:size(coexp,2)
       same=A(x,y,:);
      pval=((table2array(coexp(x,y))-mean(same))/std(same));
      if x==y
         DIFF(x,y)=10; 
      else
      DIFF(x,y)=pval;
      end
    end
end

figure
h=heatmap(coexp.Properties.VariableNames,coexp.Properties.VariableNames,DIFF);
h.Colormap=jet(30);
colormapeditor
h.Title='Z-score'





OUTPUTT=array2table(DIFF);
OUTPUTT.Properties.RowNames=coexp.Properties.VariableNames;
OUTPUTT.Properties.VariableNames=coexp.Properties.VariableNames;






end