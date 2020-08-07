% bin spatial data into hexagonal bins and output count files
% pool all data from files listed in the batch_input_file
% Xiaoyan, 2017
function[OUTPUT] =ISS_OverlappingBins_moresamples(SUPEROBJECT,distance,radius)

%% modify here
%batch_input_file = 'G:\DIPG pciseq\files_to_poolORGANOIDS.txt';
hexagon_radius =radius;   % in pixel
%output_fileprefix = 'pooled_bincounts_ORGANOIDS';
distx=distance;%90
disty=distance;%200%Number of pca components. Must be between number of samples and number of genes includeda
hexbin_size = distx*0.30;   % size of plots

samplenames=fieldnames(SUPEROBJECT);

allNames = {};
mapNames = cell(1,1);
binCounts = cell(1,1);
binPos = cell(1,1);
scales={};
images={};
for i=1:size(samplenames,1)
% start processing
images{i}=SUPEROBJECT.(samplenames{i}).image;
scales{i}=SUPEROBJECT.(samplenames{i}).scale;
name=SUPEROBJECT.(samplenames{i}).name;
pos=SUPEROBJECT.(samplenames{i}).pos;
[uNames, ~, iName] = unique(name);
[inbin, bincenters,iName] = hexbin_alldistr2(pos, hexagon_radius,distx,disty,iName);

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
end

% organize into one big matrix
rNames = [];
for i = 1:size(samplenames,1)
 rNames = [rNames,...
        repelem(i,size(binPos{i},1))];
end
%
matCounts = zeros(length(rNames), numel(allNames));
nrow = 0;
for i = 1:size(samplenames,1)
    matCounts(nrow+1:nrow+size(binPos{i}),mapNames{i}) = binCounts{i};
    nrow = nrow + size(binPos{i});    
end

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
 binPos = cat(1, binPos{:}); 
% binPos = binPos(1:end,:); %subsection
% binPosWrite = [rNames, num2cell(binPos(:,2:3))]';
% fid = fopen([output_fileprefix, '_binpos.csv'], 'w');
% fprintf(fid, 'file_bin,bincenter_x,bincenter_y\n');
% fprintf(fid, '%s,%d,%d\n', binPosWrite{:});
% fclose(fid);

mean(sum(matCounts,2));

%cor=corrcoef(matCounts);
%clustergram(cor,'ColumnLabels',uNames,'RowLabels',uNames);


OUTPUT.exp=matCounts;
OUTPUT.loc=binPos(:,[3,2]);
OUTPUT.genename=sortedNames;
OUTPUT.hexbinsize=hexbin_size;
OUTPUT.sample=rNames';
OUTPUT.image=images;
OUTPUT.scale=scales;

end