function [MeanClassExp,ClassNames,GeneNames] = expressionmat(SPOTS,gSet)
% o = o.call_cells(gSet)
%  
% Cell calling via negative binomial model
% 
% input gSet: GeneSet structure containing results of scRNAseq clustering
%
% Creates outputs: 
% pCellClass: posterior probability of each cell to be in each class (nCells*nClasses)
% pSpotCell: posterior probability of each spot to be in top 5 neighboring cells (nSpots * nCells, sparse)
% note that the last class is a zero-expressing cell; the last cell is background

o=iss;
o.Inefficiency=SPOTS.inefficiency;
%% load properties of local region of interest
%load(o.CellMapFile); % CellMap and y0, y1, x0, x1 that are its coords in full image

%% diagnostic parameters in local coordinates
% o.CellCallShowCenter = [15070 6357];
% o.CellCallShowRad = 200;
% o.ExampleCellCenter = o.CellCallShowCenter;

% exclude genes that are useless, which is none of them?
%  ExcludeGenes = {'3110035E14Rik', 'Vsnl1', 'Atp1b1', 'Slc24a2', 'Tmsb10', 'Calm2', 'Gap43', 'Fxyd6'};
Uni=unique(SPOTS.spotname);
ExcludeGenes=Uni(~ismember(Uni,gSet.GeneName));

% %LOAD
%expa=Expression(:,3:4);
expa=SPOTS.location;
%Expressiontab=readcell('I:\MIPs_ms120\sergio\CellprofilerALL\QT_0.2_details_noNNNN.csv');
%Expressiontab=readcell('K:\iss-analysis-master\lib\omero\bfmatlab\PreprocessingERIK3\Cellprofiler1\QT_0.35_details_noNNNN.csv');
expnames=SPOTS.spotname;
AllGeneNames=expnames;
o.SpotGlobalYX=expa;
% %END LOAD
o.Inefficiency=1;
%% Run spot inclusion always

IncludeSpot = ~ismember(AllGeneNames, ExcludeGenes);%...
    %& o.SpotScore> quality;
% SpotYX is only the spots we are bothered with, in global coordinates
SpotYX = o.SpotGlobalYX(IncludeSpot,:);
SpotGeneName = AllGeneNames(IncludeSpot);
SpotYX=[round(SpotYX(:,2)),round(SpotYX(:,1))];


%rp = regionprops(SPOTS.Cellmaps{1});
%CellYX = fliplr(vertcat(rp.Centroid)); % convert XY to YX
%CellArea0 = vertcat(rp.Area); 

figure(12498);
gscatter(SpotYX(:,2),SpotYX(:,1),SpotGeneName);
%hold on
%scatter(CellYX(:,2),CellYX(:,1),'r');



%% get info about cells

%MeanCellRadius = mean(sqrt(CellArea0/pi))*.5; % the dapi part is only half of the typical radius
%RelCellRadius = [sqrt(CellArea0/pi)/MeanCellRadius; 1]; % but here we want the whole thing

%% get arrays ready

% SpotGene(nS): which gene is each spot
% MeanClassExp(nK,nG): mean expression of each gene in each class
% Neighbors(nS, nN): closest neighboring cells for each spot
% D(nS, nN): distance penalty for each of these
% GeneNames(nG): name of each gene
% ClassNames(nK): name of each class

[GeneNames, ~, SpotGeneNo] = unique(SpotGeneName);
TotGeneSpots = accumarray(SpotGeneNo,1);
%gSet.Class=gSet.Class';
ClassNames = vertcat(unique(gSet.Class', 'stable'),{'Zero'});

nG = length(GeneNames);
nK = length(ClassNames); % last is zero-expression
%nC = size(CellYX,1)+1; % last is misreads
nS = size(SpotYX,1);
nN = o.nNeighbors+1; % last is misreads (always a neighbor)

ClassPrior = [.5*ones(1,nK-1)/nK .5];


ClassDisplayNames = ClassNames;

MeanClassExp = zeros(nK, nG);

%In case of adding clustering data, add here
gSet.GeneExp=gSet.GeneExp;
gSub = gSet.GeneSubset(GeneNames);
gSub.GeneExp=gSub.GeneExp';
g=gSub;
p=1;
q=1;
%THIS IS MEAN CALCULATION
for k=1:nK-1 % don't include last since it is zero-expression class
            Cells = g.IdentifyCells(ClassNames{k});
            h = g;
            h.GeneExp = g.GeneExp(Cells,:);
            h.GeneName = g.GeneName;
            h.nGenes = g.nGenes;
            h.nCells = length(Cells);
            if ~isempty(g.CellName), h.CellName = g.CellName(Cells); end
            if ~isempty(g.tSNE), h.tSNE = g.tSNE(:,Cells); end

            h.GeneExp=h.GeneExp';
            Norm = mean(h.GeneExp.^q,1).^(1/q);
            Not0Norm = (Norm>0);
            s = h.CellSubset(Not0Norm);
           % s.GeneExp=s.GeneExp';
            s.GeneExp=s.GeneExp';
            s.GeneExp = bsxfun(@rdivide, h.GeneExp(:,Not0Norm), Norm(Not0Norm).^p);
            ScaleFac = (sum(h.GeneExp(:).^q)/sum(s.GeneExp(:).^q)).^(1/q);
            s.GeneExp = s.GeneExp*ScaleFac;
    MeanClassExp(k,:) = (o.Inefficiency .* mean(s.GeneExp,2))';
end

end