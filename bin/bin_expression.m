% use bin count data to run PCA and tSNE
% Xiaoyan, 2018
function bin_expression(CELLS)


%% do not modify
% import data
tableCount = array2table(CELLS.expressionmatrix);
tableCount.Properties.VariableNames=CELLS.expressionnames;
%Filter
%r for 1 column
vect=sum(table2array(tableCount(:,2:end)),2) >= CELLS.count_threshold;
tableNO=tableCount(sum(table2array(tableCount(:,2:end)),2) < CELLS.count_threshold, :);
tableCount=tableCount(sum(table2array(tableCount(:,2:end)),2) >= CELLS.count_threshold, :);


GOOD=ismember(tableCount.Properties.VariableNames,CELLS.Genes_of_interest);


% original column names
meas=tableCount(:,GOOD);
meas=mean(table2array(meas),2);

% get position
pos = CELLS.location;
pos(:,3)=vect;
pos=array2table(pos);
posno=pos(pos.pos3==0,:);
posno=table2array(posno(:,1:2));
pos=pos(pos.pos3==1,:);
pos=table2array(pos(:,1:2));

% visualize tSNE in RGB (no background)
%pos=vertcat(pos,posno);
tab=zeros(size(posno,1),3);
tab(:,:)=1;


samplex=unique(CELLS.sample);
SUBSAMPLE=CELLS.sample(vect);
SUBCOORD=CELLS.location(vect,:);

for element=1:size(samplex,1)
elename=samplex(element);
DAP=imread(char(CELLS.images(element)));
DAP=imresize(DAP,1/cell2mat(CELLS.scale(element)));
if size(DAP,3)>1
    DAP=rgb2gray(DAP);
end
 Y=[meas(ismember(SUBSAMPLE,elename),:), meas(ismember(SUBSAMPLE,elename),:)./meas(ismember(SUBSAMPLE,elename),:), meas(ismember(SUBSAMPLE,elename),:)];
 SEC=rgbscale(Y);
figure;
imshow(DAP*5);
hold on
s=scatter(SUBCOORD(ismember(SUBSAMPLE,elename),1),SUBCOORD(ismember(SUBSAMPLE,elename),2),CELLS.hexbin_size,SEC,'filled','Marker', 's');
s.MarkerFaceAlpha = 0.70;
hold off
end




end