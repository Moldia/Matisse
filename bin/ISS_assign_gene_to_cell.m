function [TOTAL,EXPRESSIONMAT,CELS]= ISS_assign_gene_to_cell(SPOTS,CellMap,XY)
XY=fliplr(vertcat(regionprops(CellMap).Centroid));
details4sel(:,2)=SPOTS.name;
details4sel(:,3)=num2cell(SPOTS.pos(:,1));
details4sel(:,4)=num2cell(SPOTS.pos(:,2));

VAL=[];
DIS=[];
for sp = 1:size(details4sel,1)
 M=sqrt(((cell2mat(details4sel(sp,4))-XY(:,1)).^2)+((cell2mat(details4sel(sp,3))-XY(:,2)).^2));
 VAL=[VAL,find(M==min(M))];
 DIS=[DIS,min(M)];
end

figure
scatter((cell2mat(details4sel(:,4))),(cell2mat(details4sel(:,3))));
hold on 
scatter(XY(:,1),XY(:,2));

%make all variables together in a table
TOTAL={};
TOTAL(:,1)=num2cell(DIS');
TOTAL(:,2)=num2cell(VAL');
TOTAL(:,3)=details4sel(:,2);
TOTAL(:,4)=details4sel(:,3);
TOTAL(:,5)=details4sel(:,4);


CELS=unique(cell2mat(TOTAL(:,2)));
EXPRESSIONMAT.exp=crosstab(cell2mat(TOTAL(:,2)),(TOTAL(:,3)));
EXPRESSIONMAT.loc=XY(CELS,:);
EXPRESSIONMAT.hexbinsize=10;
EXPRESSIONMAT.genename=unique(SPOTS.name);
end