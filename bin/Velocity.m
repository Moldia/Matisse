% plotting on top of the density estimate
% Xiaoyan, 2020
function SPOTS = Velocity(SPOTS)

%% modify here
sam=unique(SPOTS.sample);

for eachsam=1:size(sam,1)

sampl=sam(eachsam);
SELECTED=ismember(SPOTS.sample,sampl);
name=SPOTS.spotname(SELECTED,:);
pos=SPOTS.location(SELECTED,:);
image=SPOTS.images{eachsam};
scale=SPOTS.scale{eachsam};
%gene_density=SPOTS.Genes_of_interest{1};
bandwid=SPOTS.bandwidth;
% kde

SPO=unique(SPOTS.spotname);




for element=1:size(SPO,1)
disp(element);
gene_density=SPO(element);
bandwid=SPOTS.bandwidth;
density = gene_kde(name, pos, gene_density, bandwid, image, scale);
XSPOTS=[(1+SPOTS.distance):SPOTS.distance:(size(density,2)-SPOTS.distance)];
YSPOTS=[(1+SPOTS.distance):SPOTS.distance:(size(density,1)-SPOTS.distance)];

VSPOTS=[];
USPOTS=[];
XREAL=[];
YREAL=[];
for XS=XSPOTS
for YS=YSPOTS
    XREAL=[XREAL,XS];
    YREAL=[YREAL,YS];
    VSPOTS=[VSPOTS,(density(YS-SPOTS.distance,XS)- density(YS+SPOTS.distance,XS))];
    USPOTS=[USPOTS,(density(YS,XS-SPOTS.distance)- density(YS,XS+SPOTS.distance))];
end
end
% figure
% quiver(XREAL,YREAL,USPOTS,VSPOTS)

if element==1
   MATRIXR=zeros(size(XREAL,2),4,size(SPO,1));
end

case1=(VSPOTS<0 & USPOTS<0);
VSPOTS(case1)=-VSPOTS(case1);
USPOTS(case1)=-USPOTS(case1);

case2=(VSPOTS>0 & USPOTS<0);
VSPOTS(case2)=-VSPOTS(case2);
USPOTS(case2)=-USPOTS(case2);

MATRIXR(:,1,element)=abs(XREAL);
MATRIXR(:,2,element)=abs(YREAL);
MATRIXR(:,3,element)=VSPOTS;
MATRIXR(:,4,element)=USPOTS;


% % plot
% figure;
% imshow(density,[]);
% hold on;
% colormap(gca, SPOTS.colormap);
% title([gene_density ' density']);
end


%HERE WE ADD FINAL CALCS
MATFIN=sum(MATRIXR,3);
figure
imshow(imresize(imread(SPOTS.images{eachsam}),0.2)*4);
hold on 
quiver(MATRIXR(:,1,1),MATRIXR(:,2,1),MATFIN(:,3)*100,MATFIN(:,4)*100)

end
end 
