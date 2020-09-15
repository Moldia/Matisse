%% iss: code for processing of in situ sequencing
% Designed by Sergio Marco Salas
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
%
% to use:
% spots = matisse; % create structure, default parameters
% % change any parameters you want, and set o.FileN

% this current file MATISSE.m just contains default parameter values

classdef matisseMOD
    properties
        %% interactive graphics mode. 1 means some, 2 means a lot.
        Graphics = 1;
        
        %% Create empty variables that have to be filled
        %Spot name includes the name of each spot/cell
        spotname=[]; 
        %Location includes the location in XY of each spot/cell
        location=[];
        %Sample includes the sample where each spot/cell belongs to
        sample=[];
        %Domain includes the domain that one cell has been assigned to
        domain=[];
        %he unique samples that we have
        sample_ID=[];
        %Images includes the background images that we have
        images=[];
        %Scale includes the scale of each image
        scale=[];
        %% Expression variables
        %If we are working with cells, expression matrix includes the
        %expression of cell
        expressionmatrix=[];
        %In case of working with expression matrix, this includes the genes
        %on each row
        expressionnames=[];
        centroidpos={};
        Cellmaps={};
        centroidID={};
       
        %% Parameters
        %This parameters are function-specific and can be used as default o
        % modify them
        
        %GeneDensity parameters
       % GeneDensity;
        GeneDensity_Genes_of_interest={};
        GeneDensity_bandwidth=300;
        GeneDensity_colormap=parula;
        %OverlayTwoDensities
        OverlayTwoDensities_Genes_of_interest={};
        OverlayTwoDensities_bandwidth=300;
        OverlayTwoDensities_colormap=parula;
        
        %COLOCALIZATION
        Colocalization_number_neighbors=3;
        
        % ROI
        ROI_bandwidth=300;
        
        % Gradients parameters
        Gradients_limit=10000;
        Gradients_bandwidth=200;
        Gradients_colormap=parula;
        
        %Find_gradients
        Find_gradients_bandwidth=200;
        Find_gradients_distance=20;
        
        %Bins
        Bins_distance=5000;
        Bins_radius=5000;
        
        %Principal_components
        Principal_components_colormap=redbluecmap;
        
        %Filter_cells
        FC_MIN_counts_cell=5;
        FC_MAX_counts_cell=200;
        FC_MIN_counts_gene=20;
        
        %Normalize_cells
        NC_method='log2'; % log2 / standardize
        
        %CLUSTERING
        SpotEnvBin_radius=100;
       
        %DOMAIN
        
        
        
        %%%%
        bandwidth=200;
        hexbin_size=300;
        radius=300;
        distance=300; %this is needed for Binning
        Genes_of_interest={};
        colormap=parula;
        count_threshold=5;
        %Segmentation parameters
        segmentation_percent=95;
        
        %IDs for further processes
        Clustering_ID=[];
        RGB_CODE=[];
        UMAP_clustering=[];
        pciseq_clustering=[];

        %Gradients
        limit=800;
        
        %PCISEQ
        pciseq_probabilities=[];
        pciseq_probabilities_names=[];
        CellCallRegionYX=[];
        inefficiency=0.02;
        
    end
    
end

 





