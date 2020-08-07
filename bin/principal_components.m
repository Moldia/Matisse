function[NEWMAT,correlat]= principal_components(CELLS,varargin)

if nargin==1
imagesel=1;
else
 imagesel=varargin{1}
end


allsamples=unique(CELLS.sample);
SELECTED=ismember(CELLS.sample,CELLS.sample);

SELECTIONS=CELLS;
SELECTIONS.spotname=SELECTIONS.spotname(SELECTED,:);
SELECTIONS.location=SELECTIONS.location(SELECTED,:);
SELECTIONS.sample=SELECTIONS.sample(SELECTED,:);
SELECTIONS.images=SELECTIONS.images;
SELECTIONS.scale=SELECTIONS.scale(imagesel);
SELECTIONS.expressionmatrix=SELECTIONS.expressionmatrix(SELECTED,:);



[coeff,score,latent,~,explained] = pca(SELECTIONS.expressionmatrix);
figure
bar(explained)
ylabel('Percentage of variable explained');
xlabel('Number of components');
title('Explained variability by different PC')
correlat=zeros(10,size(SELECTIONS.expressionmatrix,2));
for pc=1:10
    for var=1:size(SELECTIONS.expressionmatrix,2)
       correlat(pc,var)= corr(score(:,pc),SELECTIONS.expressionmatrix(:,var));
    end
end

correlat=array2table(correlat');
correlat.Properties.RowNames=SELECTIONS.expressionnames;


figure;
heatmap(table2array(correlat),'Colormap',CELLS.Principal_components_colormap,'GridVisible','off',...
    'YDisplayLabels',correlat.Properties.RowNames); % 'XDisplayLabels',correlat.Properties.RowNames,'YDisplayLabels',DE

allsamples=unique(CELLS.sample);


LOCATIONS=[];
for vis=1:size(allsamples,1);
SELECTED2=ismember(CELLS.sample,allsamples(vis));
SELECTIONS2=CELLS;
SELECTIONS2.spotname=SELECTIONS.spotname(SELECTED2,:);
SELECTIONS2.location=SELECTIONS.location(SELECTED2,:);
SELECTIONS2.sample=SELECTIONS.sample(SELECTED2,:);
SELECTIONS2.images=SELECTIONS.images;
SELECTIONS2.scale=SELECTIONS.scale(imagesel);
SELECTIONS2.expressionmatrix=SELECTIONS.expressionmatrix(SELECTED2,:);
SECE=score(SELECTED2,:);

if size(LOCATIONS,1)>0
LOCATIONS=[LOCATIONS;SELECTIONS2.location(:,1)+max(LOCATIONS(:,2))+2000,(SELECTIONS2.location(:,2))];
SECES=[SECES;SECE];
else
LOCATIONS=[SELECTIONS2.location(:,1),(SELECTIONS2.location(:,2))];
SECES=[SECE];
end

figure
%I=imresize(imread(CELLS.images{imagesel}),1/CELLS.scale{imagesel});
for CAL=1:10
subplot(2,5,CAL);
VE=SECES(:,CAL);
a=scatter(LOCATIONS(:,1),LOCATIONS(:,2),round((20000/size(LOCATIONS,1)))+5,round(VE/std(VE)),'square','filled');
hold on
%h=imshow(I*20);
%alpha(0.05);
colormap(CELLS.Principal_components_colormap);
title([,'PC',num2str(CAL)]);
linkaxes()
%set(a,'Color','k')
end

end

DEVS=[];
if size(unique(CELLS.sample),1)>1
 for pc=1:size(score,2)
    disp(pc);
    [B,dev,stats] = mnrfit(score(:,pc),categorical(CELLS.sample));
    DEVS=[DEVS,dev];
 end 
end

figure
bar(max(DEVS)-DEVS,'r');
ylabel('Variance explained by sample');
xlabel('Number of components');
title('Explained variability by different PC')

prompt = 'List the selected Principal components';
x = input(prompt)

NEWMAT=SELECTIONS;
NEWMAT.expressionmatrix=score(:,[x])*coeff([x],:);

end