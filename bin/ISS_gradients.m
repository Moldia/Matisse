
function [] = ISS_gradients(SPOTS,limit,bandwidth)
%This script is the first part in gradient analysis for circular tissues
%(from outside to inside)

%This specify both the basis image and the details of the gene
I=imread(SPOTS.image);
%decoded_file='L:\Organoid_method_development\Organoids_erik\Organoid3_Output\7\QT_0.6_details_noNNNN.csv';
scale=1; %SALE OF THE DAPI
roi_number=1; %NUMBER OF ROIs
%limit=800; % Maximum distance to measure the gradient
%bandwidth=30; %Bandidth 


%%%%%%%% DO NOT MODIFY FROM HERE
%table=readtable(decoded_file);

[Coord, Coord_write, blank] = drawpolygon(I, roi_number, scale);


%DENSITY OUT
dist=[];
for qul =1:size(SPOTS.pos,1)
disp(qul)
sele=SPOTS.pos(qul,:);
MAT=sqrt(min((Coord_write(:,2)-sele(:,1)).^2+(Coord_write(:,3)-sele(:,2)).^2));
dist=[dist,MAT]; 
    
end

distance={};
distance=cell2table(distance);
distance.Gene=SPOTS.name;
distance(:,2:3)=num2cell(SPOTS.pos);
distance.Density=dist';

%Plot 
figure
imshow(imresize(imread(SPOTS.image),1/SPOTS.scale)*1);
hold on 
scatter(table2array(distance(:,2)),table2array(distance(:,3)),[],distance.Density*10);
colormap(hot);

cat=unique(distance.Gene);
figure;
TAC=[];
DIS=[];
NAM=[];
for gens=1:size(cat);
gen=cat(gens);
distsel=distance(ismember(distance.Gene,gen),:);
as=fitdist(distsel.Density,'kernel','BandWidth',bandwidth);
set=pdf(as,1:1:limit);
si=plot(1:1:limit,set);
DIS=[DIS,si];
si;
hold on
end
legend(DIS,cat);

