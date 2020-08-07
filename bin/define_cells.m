%Function: define_cells
%2020, Nilsson Lab, Scilifelab
%Objective: Combine the cell position with the gene position to define
%cells and cell expression, saved in a NEW matisse object containing cells
%Input: Matisse object, containing a cell map (after DAPI_segmentation)
%Output: Matisse object, that now contains the Cells instead of spots

function [CELLS]= define_cells(SPOTS)

%We first create a CELLS object, based on matisse.m
CELLS=matisseMOD;
sam=unique(SPOTS.sample);

XYLOC=[];
CENAME=[];
GENAMES=[];
MATRIX=[];
SAMPLE=[];

for eachsam=1:size(sam,1)
sampl=sam(eachsam);
SELECTED=ismember(SPOTS.sample,sampl);
names=SPOTS.spotname(SELECTED,:);
pos=SPOTS.location(SELECTED,:);


XY=SPOTS.centroidpos{eachsam};
details4sel={};
details4sel(:,2)=names;
details4sel(:,3)=num2cell(pos(:,1));
details4sel(:,4)=num2cell(pos(:,2));
%details4sel=details4sel(SELECTED,:);
VAL=[];
DIS=[];
disp(size(details4sel));
for sp =  1:size(details4sel,1)
 M=sqrt(((cell2mat(details4sel(sp,4))-XY(:,1)).^2)+((cell2mat(details4sel(sp,3))-XY(:,2)).^2));
ce=find(M==min(M));
 VAL=[VAL,ce(1)];
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
[EXPRESSIONMAT.exp,chi2,p,labels]=crosstab(cell2mat(TOTAL(:,2)),(TOTAL(:,3)));
EXPRESSIONMAT.loc=XY(CELS,:);
EXPRESSIONMAT.hexbinsize=10;
EXPRESSIONMAT.genename=labels(1:size(EXPRESSIONMAT.exp,2),2);
mate=array2table(EXPRESSIONMAT.exp);
mate.Properties.VariableNames=EXPRESSIONMAT.genename;
NONUSED=unique(SPOTS.spotname(~ismember(SPOTS.spotname,mate.Properties.VariableNames)));
for s=1:size(NONUSED)
   vTbl = table(zeros(size(mate,1),1), 'VariableNames',NONUSED(s)); 
   mate=[mate,vTbl]; 
   mate=mate(:,sort(mate.Properties.VariableNames));
end


XYLOC=[XYLOC;EXPRESSIONMAT.loc];
CENAME=[CENAME;CELS];
SAMPLE=[SAMPLE,repelem(sampl,size(EXPRESSIONMAT.loc,1))];
GENAMES=labels(1:size(EXPRESSIONMAT.exp,2),2);
MATRIX=[MATRIX;table2array(mate)];

end

CELLS.location=[XYLOC(:,2),XYLOC(:,1)];
CELLS.spotname=CENAME;
CELLS.expressionmatrix=MATRIX;
CELLS.expressionnames=mate.Properties.VariableNames';
CELLS.sample=SAMPLE';

CELLS.images=SPOTS.images;
CELLS.scale=SPOTS.scale;

end