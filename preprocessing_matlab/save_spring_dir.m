function save_spring_dir(E,D,k,gene_list,project_directory, varargin)
% Create a SPRING project directory and write files to be loaded by SPRING
% 
% Required inputs
%   E                  = Matrix of gene expression. Rows correspond to
%                        cells and columns correspond to genes. 
%   D                  = Distance matrix for construction of knn graph.
%                        Any distance matrix can be used as long as higher
%                        values correspond to greater distances.
%   k                  = Number of edges assigned to each node in knn graph
%   gene_list          = An ordered cell array of gene names. The cell 
%                        array should have length size(E,2)
%   project_directory  = Path to a directory where SPRING readable files
%                        will be written. The directory does not have to 
%                        exist before running this function.
% Optional inputs
%   "custom_colors"    = cell array with one row for each custom color 
%                        track. The first entry in each row is the color 
%                        tracks and N cells, this should be a T x (N+1) 
%                        cell array. 
%   "cell groupings"   = cell array with one row for each cell grouping. 
%                        The first entry in each row is the name of the
%                        of the grouping (e.g. "sampleID") and each 
%                        subsequet entry is a cell label (e.g. "sample_1").
%                        If there are T different groupungs and N cells, 
%                        then this should be a T x (N+1) cell array. 
%
%%
    cell_groupings = cell(0);
    custom_colors = cell(0);
    vargs = varargin;
    nargs = length(vargs);
    names = vargs(1:2:nargs);
    values = vargs(2:2:nargs);
    for i=1:length(names)
        if strcmp(names{i},'cell_groupings')
            cell_groupings = values{i};
        end 
        if strcmp(names{i},'custom_colors')
            custom_colors = values{i};
        end
    end
    for i=1:size(custom_colors,1)
        for j=2:size(custom_colors,2)
            val = custom_colors{i,j};
            if ischar(custom_colors{i,j})
                custom_colors{i,j} = str2double(custom_colors{i,j});
            end
        end
    end
               
    % Prepare output directory
    mkdir(project_directory)
    if project_directory(end) ~= '/'
        project_directory = strcat(project_directory,'/');
    end
    
    % Build graph
    disp('Building graph');
    edges = get_knn_edges(D,k);
    write_graph(size(D,1),edges,project_directory);
    
    % Save custom color tracks
    disp('Writing custom colors');
    write_color_tracks(custom_colors,strcat(project_directory,'color_data_gene_sets.csv'));
    
    % Save gene colortracks
    disp('Writing gene colors');
    mkdir(strcat(project_directory,'gene_colors'));
    II = round(length(gene_list) / 50) + 1;
    for j=1:50
        disp(strcat('Writing block_',int2str(j)));
        fname = strcat(project_directory,'gene_colors/color_data_all_genes-',num2str(j-1),'.csv');
        all_gene_colors = cell(1,size(E,1)+1);
        current_row = 1;
        for i = (II*j+1):(II*(j+1)+1)
            if i <= size(E,2) && mean(E(:,i)) > 0.05       
                all_gene_colors{current_row,1} = gene_list{i};
                for k = 1:size(E,1)
                    all_gene_colors{current_row,k+1} = E(k,i);
                end
                current_row = current_row + 1;
            end
        end
        write_color_tracks(all_gene_colors, fname);
    end
    
    % Save color tracks summary statistics
    disp('Saving color track statistics');
    colorstats = cell(1);
    fieldlist = cell(1);
    for i = 1:size(E,2)
        stats = [mean(E(:,i)), std(E(:,i)), 0, max(E(:,i))+.01, prctile(E(:,i),99)];
        colorstats{i} = stats;
        fieldlist{i} = gene_list{i};
    end
    for i = 1:size(custom_colors,1)
        track = cell2mat(custom_colors(i,2:end));
        stats = [mean(track),std(track),0, max(track)+.01, prctile(track,99)];
        colorstats{i+size(E,2)} = stats;
        fieldlist{i+size(E,2)} = custom_colors{i,1};
    end
    colorstats = cell2struct(colorstats,fieldlist,2);  
    json.write(colorstats, strcat(project_directory,'color_stats.json'));
    
	% Save categorical coloring data
    disp('Saving categorical coloring data');
    categorical_coloring_data = struct();
    for i=1:size(cell_groupings,1)
        labelset = unique(cell_groupings(i,2:end));
        label_colors = cell2struct(get_categorical_colors(length(labelset)),labelset,1);
        categorical_coloring_data.(cell_groupings{i,1}) = struct('label_colors',label_colors,'label_list',char(cell_groupings(i,2:end)));
    end
    json.write(categorical_coloring_data,strcat(strcat(project_directory,'categorical_coloring_data.json')));
end