%Function: OverlappingBins
%2020, Nilsson Lab, Scilifelab
%Objective: Bin the spots distributed in the tissue in several bins
%distributed at certain distances for further analysis of this bins
%Input: Matisse object 
  %NOTE THAT: parameters such as SPOTS.Bins_radius or SPOTS.Bins_distance can be
  %previusly modified in the matisse object
%Output: Matisse objects reprenting each bin. 
function[CELLS] =OverlappingBins(SPOTS)

hexagon_radius =SPOTS.Bins_radius;   % in pixel
distx=SPOTS.Bins_distance;%90
disty=SPOTS.Bins_distance;%200%Number of pca components. Must be between number of samples and number of genes includeda
hexbin_size = distx*0.40;   % size of plots
    
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
name=SPOTS.spotname(SELECTED,:);
pos=SPOTS.location(SELECTED,:);

allNames = {};
mapNames = cell(1,1);
binCounts = cell(1,1);
binPos = cell(1,1);


    [uNames, ~, iName] = unique(name);
    [inbin, bincenters,iName] = hexbin_alldistr2(pos, hexagon_radius,distx,disty,iName);
    i=1
    % count transcripts in each bin
    counts = histcounts2(inbin, iName ,...
        [unique(inbin(~ismember(inbin,0)))', max(inbin)+1],...
        1:numel(uNames)+1);
    binCounts{i} = counts;
    binPos{i} = [unique(inbin(~ismember(inbin,0))), bincenters];
    
    % add genes that have not appeared yet
    allNames = [allNames; setdiff(uNames, allNames)];
        
    % name index in the pool
    iuNames = cellfun(@(v) find(strcmp(v, allNames)), uNames);
    mapNames{i} = iuNames;    


% organize into one big matrix
rNames = {};
rNames = [rNames;catstrnum([strtok(SPOTS.images{eachsam}, '.'), '_r', num2str(hexagon_radius), '_hexbin'], binPos{i}(:,1))];

%
matCounts = zeros(length(rNames), numel(allNames));
nrow = 0;
    matCounts(nrow+1:nrow+size(binPos{i}),mapNames{i}) = binCounts{i};
    nrow = nrow + size(binPos{i});
[sortedNames, orderNames] = sort(allNames);
matCounts = matCounts(:,orderNames);

rNames=rNames(1:end); %subsection
matCounts=matCounts(1:end,:); %subsection

% % write count file 
% binCountsWrite = [rNames, num2cell(matCounts)]';
% fid = fopen([output_fileprefix, '_count.csv'], 'w');
% fprintf(fid, 'file_bin,');
% fprintf(fid, lineformat('%s', numel(sortedNames)), sortedNames{:});
% fprintf(fid, ['%s,' lineformat('%d', numel(sortedNames))], binCountsWrite{:});
% fclose(fid);


% write hexbin position file
% binPos = cat(1, binPos{:}); 
% binPos = binPos(1:end,:); %subsection
% binPosWrite = [rNames, num2cell(binPos(:,2:3))]';
% fid = fopen([output_fileprefix, '_binpos.csv'], 'w');
% fprintf(fid, 'file_bin,bincenter_x,bincenter_y\n');
% fprintf(fid, '%s,%d,%d\n', binPosWrite{:});
% fclose(fid);
%cor=corrcoef(matCounts);
%clustergram(cor,'ColumnLabels',uNames,'RowLabels',uNames);

EXPRESSIONMAT.exp=matCounts;
binPos=binPos{1,1};
EXPRESSIONMAT.loc=binPos(:,[3,2]);
EXPRESSIONMAT.genename=uNames;
EXPRESSIONMAT.hexbinsize=hexbin_size;

mate=array2table(EXPRESSIONMAT.exp);
for ex= 1:size(EXPRESSIONMAT.genename )
    if ismember(EXPRESSIONMAT.genename(ex),'')
     EXPRESSIONMAT.genename(ex)={'empty'}
        
        
    end
end
mate.Properties.VariableNames=EXPRESSIONMAT.genename;

for ex= 1:size(SPOTS.spotname,1 );
    if ismember(SPOTS.spotname(ex,1),'');
     SPOTS.spotname(ex,1)={'empty'}  ; 
    end
end
NONUSED=unique(SPOTS.spotname(~ismember(SPOTS.spotname,mate.Properties.VariableNames)));


for s=1:size(NONUSED)
   vTbl = table(zeros(size(mate,1),1), 'VariableNames',NONUSED(s)); 
   mate=[mate,vTbl]; 
   mate=mate(:,sort(mate.Properties.VariableNames));
end

XYLOC=[XYLOC;EXPRESSIONMAT.loc];
CENAME=[CENAME;rNames];
SAMPLE=[SAMPLE,repelem(sampl,size(EXPRESSIONMAT.loc,1))];
GENAMES=unique(SPOTS.spotname);
MATRIX=[MATRIX;table2array(mate)];
end


CELLS.location=[XYLOC(:,2),XYLOC(:,1)];
CELLS.spotname=CENAME;
CELLS.expressionmatrix=MATRIX;
CELLS.expressionnames=GENAMES;
CELLS.sample=SAMPLE';

CELLS.images=SPOTS.images;
CELLS.scale=SPOTS.scale;
CELLS.hexbin_size=hexbin_size;
end