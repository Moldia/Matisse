function[]=ISS_TSNE_heterogeneity(EXPRESSIONMAT,SPOTS,minimum_expression,hexbin_size)
% import data
tableCount=array2table(EXPRESSIONMAT.exp);

%Filter for 1 column
vect=sum(table2array(tableCount(:,2:end)),2) >= minimum_expression;
tableNO=tableCount(sum(table2array(tableCount(:,2:end)),2) < minimum_expression, :);
tableCount=tableCount(sum(table2array(tableCount(:,2:end)),2) >= minimum_expression, :);





% original column names
cNames = EXPRESSIONMAT.genename;

% make sure there are no multiple entries of the same gene
assert(numel(cNames)== numel(unique(cNames)),...
    'Column names are not unique!')

% create checkboxes and get selected values
cbValues = checkboxes(cNames);
idx = cellfun(@(v) find(strcmp(v, cNames)), cbValues(:,1));
isSelected = false(numel(idx), 1);
isSelected(idx) = cell2mat(cbValues(:,2));
cGenes = table2array(tableCount(:,isSelected));
genes = cNames(isSelected)'

% PCA
[coeff, score, latent] = pca(cGenes);
% visualize first two components
figure, biplot(coeff(:,1:2), 'Scores', score(:,1:2), 'VarLabels', genes);
title('top two principle components');



%Filter for 1 variable



% tSNE in MATLAB
% ONLY >=R2018a
seeds = 1e-4*randn(size(cGenes,1), 3);
Y = tsne(cGenes, 'NumDimensions', 3, 'NumPCAComponents', 4, 'Perplexity', 40,...
    'Standardize', 1, 'LearnRate', 1000, 'Verbose', 1, 'InitialY', seeds); 



% get position
pos=EXPRESSIONMAT.loc;
pos(:,3)=vect;
pos=array2table(pos);
posno=pos(pos.pos3==0,:);
posno=table2array(posno(:,1:2));
pos=pos(pos.pos3==1,:);
pos=table2array(pos(:,1:2));
% visualize tSNE in RGB (no background)
%if ~hexbin_size;	hexbin_size = 10;    end
Yrgb = rgbscale(Y);

%Visualization depending on tsne dimension
figure, scatter3(Y(:,1),Y(:,2),Y(:,3),10, Yrgb,'filled');
colormap(jet(20));
title({'tSNE dim reduction to three', 'Color based on axis, according to spatial representation'});



% Visualisation depending on gene expression

for s=2:size(tableCount,2)
genc=table2array(tableCount(:,s));
ti=tableCount.Properties.VariableNames(s);
figure, scatter3(Y(:,1),Y(:,2),Y(:,3),10,genc);
colormap(hsv(100));
title(ti);
end 

%pos=vertcat(pos,posno);
tab=zeros(size(posno,1),3);
tab(:,:)=1;
%Yrgb=vertcat(Yrgb,tab)
% 
% for s = 1:numel(uSamples)
%     drawnow
%     figure(1212);
%     for elem=2000:3000%size(pos,1)
%     disp(elem);
%     hold on 
%     scatter(pos(elem,1),pos(elem,2),30, Yrgb(elem,:),'filled');
%     end
% end
% 

DAP=imread(SPOTS.image);
DAP=imresize(DAP,1/SPOTS.scale);

figure(96336);
imshow(DAP*20);
hold on
s=scatter(pos(:,2),pos(:,1),hexbin_size, Yrgb(:,:),'filled','Marker', 's');
s.MarkerFaceAlpha = 0.7;

%%
%cellclass=readcell('G:/DIPG pciseq/details_cells_DIPG_No_NNNN.csv');
%gscatter(cell2mat(cellclass(2:end,3)),cell2mat(cellclass(2:end,4)),cellclass(2:end,2),[],[],15);
% %CLUSTERING
% eucD = pdist([Yrgb],'euclidean');
% clustTreeEuc = linkage(eucD,'average'); 
% cophenet(clustTreeEuc,eucD)
% 
% figure(12121);
% [h,nodes] = dendrogram(clustTreeEuc,0);
% h_gca = gca;
% h_gca.TickDir = 'out';
% h_gca.TickLength = [.002 0];
% h_gca.XTickLabel = []; 
% hidx = cluster(clustTreeEuc,'criterion','distance','cutoff',500);
% pointsize=120;
% figure(11111);
% scatter(pos(:,1),pos(:,2),pointsize,hidx,'filled','Marker', 's')
% colormap(hsv)   
% 
% 


end