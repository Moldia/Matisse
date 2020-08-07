% bin spatial data into hexagonal bins and output count files
% pool all data from files listed in the batch_input_file
% Xiaoyan, 2017
function[OUTPUT] =ISS_OverlappingBins(SPOTS,distance,radius)

%% modify here
%batch_input_file = 'G:\DIPG pciseq\files_to_poolORGANOIDS.txt';
hexagon_radius =radius;   % in pixel
%output_fileprefix = 'pooled_bincounts_ORGANOIDS';
distx=distance;%90
disty=distance;%200%Number of pca components. Must be between number of samples and number of genes includeda
hexbin_size = distx*0.30;   % size of plots


% start processing
allNames = {};
mapNames = cell(1,1);
binCounts = cell(1,1);
binPos = cell(1,1);
name=SPOTS.name;
pos=SPOTS.pos;

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
rNames = [rNames;catstrnum([strtok(SPOTS.image, '.'), '_r', num2str(hexagon_radius), '_hexbin'], binPos{i}(:,1))];

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

mean(sum(matCounts,2));

%cor=corrcoef(matCounts);
%clustergram(cor,'ColumnLabels',uNames,'RowLabels',uNames);


OUTPUT.exp=matCounts;
binPos2=binPos{1};
OUTPUT.loc=binPos2(:,[3,2]);
OUTPUT.genename=uNames;
OUTPUT.hexbinsize=hexbin_size;

end