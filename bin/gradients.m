%Function: gradients
%2020, Nilsson Lab, Scilifelab
%Objective: Represent the gradients or 1D KDE of genes/cells in a tissue
%Input: Matisse object + sample selected (eg. 1)  
  %NOTE THAT: starting points of the gradients will need to be set by
  %selecting the limits in a ROI
%Output: Image representing gradient with X axis=distance from spots
%selected
function [varargout] = gradients(SPOTS,varargin)
if nargin==1
imagesel=1;
else
 imagesel=varargin{1}
end

roi_number=1;
allsamples=unique(SPOTS.sample);
SELECTED=ismember(SPOTS.sample,allsamples(imagesel));
SELECTIONS=SPOTS;
SELECTIONS.spotname=SELECTIONS.spotname(SELECTED,:);
SELECTIONS.location=SELECTIONS.location(SELECTED,:);
SELECTIONS.sample=SELECTIONS.sample(SELECTED,:);
SELECTIONS.images=SELECTIONS.images;
SELECTIONS.scale=SELECTIONS.scale(imagesel);

name=SELECTIONS.spotname;
pos=SELECTIONS.location;
background_image=SELECTIONS.images{1};
scale=SELECTIONS.scale{1};
%This specify both the basis image and the details of the gene
I=imread(background_image);
%decoded_file='L:\Organoid_method_development\Organoids_erik\Organoid3_Output\7\QT_0.6_details_noNNNN.csv';
roi_number=1; %NUMBER OF ROIs
limit=SPOTS.Gradients_limit; % Maximum distance to measure the gradient
%bandwidth=30; %Bandidth 


%%%%%%%% DO NOT MODIFY FROM HERE
%table=readtable(decoded_file);

[Coord, Coord_write, blank] = drawpolygon(I, roi_number, scale);


%DENSITY OUT
dist=[];
for qul =1:size(pos,1)
disp(qul)
sele=pos(qul,:);
MAT=sqrt(min((Coord_write(:,2)-sele(:,1)).^2+(Coord_write(:,3)-sele(:,2)).^2));
dist=[dist,MAT]; 
    
end

distance={};
distance=cell2table(distance);
distance.Gene=name;
distance(:,2:3)=num2cell(pos);
distance.Density=dist';

%Plot 
figure
imshow(imresize(imread(background_image),1/scale)*1);
hold on 
scatter(table2array(distance(:,2)),table2array(distance(:,3)),[],distance.Density*10);
colormap(SPOTS.Gradients_colormap);

cat=unique(distance.Gene);
RESULTS=zeros(limit,size(cat,1));
figure;
TAC=[];
DIS=[];
NAM=[];
figure
for gens=1:size(cat);
gen=cat(gens);
title('Gene density distribution along the gradient');
distsel=distance(ismember(distance.Gene,gen),:);
as=fitdist(distsel.Density,'kernel','BandWidth',SPOTS.Gradients_bandwidth);
sete=pdf(as,1:1:limit);
sete=(sete./max(sete));
si=plot(1:1:limit,sete);
%si=si./max(si); %This line has been added
DIS=[DIS,si];
si;
hold on
RESULTS(:,gens)=sete;
end
if strcmp(class(cat),'double')
    cat=num2str(cat);
   legend(DIS',cat);
else
legend(DIS,cat);
end

%%%%%%%%THIS IS FOR GRADIENT WITH SPOTS SPLITTED IN PLOTS%%%%%%%%%%%
RESULTS=array2table(RESULTS);
RESULTS.Properties.VariableNames=cat;
varargout{1}=RESULTS;
figure
title('Spot distribution along the gradient');
for gens=1:size(cat);
gen=cat(gens);
distsel=distance(ismember(distance.Gene,gen),:);
%as=fitdist(distsel.Density,'kernel','BandWidth',SPOTS.Gradients_bandwidth);
%set=pdf(as,1:1:limit);
%set=(set./max(set));
subplot(round(size(cat,1)/4)+1,4,gens)
scatter(distsel.Density,rand(size(distsel.Density,1),1).*100,6,'r','filled');
xlim([0,max(distance.Density)]);
ylabel(gen);
set(gca,'YTick',[]);
set(get(gca,'YLabel'),'Rotation',0,'VerticalAlignment','middle');
if gens~=size(cat,1)
  set(gca,'XTick',[]);
end    
%si=si./max(si); %This line has been added
end

end
