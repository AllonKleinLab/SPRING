function save_spring_dir(E,D,k,gene_list,project_directory, varargin)
% Create a SPRING project directory and write files to be loaded by SPRING
% 
% write_color_tracks(ctracks, fname)
%
% ctracks  = Cell array. Each row is a color track. The first entry in the
%            row is the name of the track and subsequent entries give 
%            specific color values. If there are N cells, then this cell
%            array should have N+1 columns. 
% fname    = Path to save the color tracks file. 
%
%%
    cell_groupings = cell(0);
    custom_colors = cell(0);
    vargs = varargin;
    nargs = length(vargs);
    names = vargs(1:2:nargs);
    values = vargs(2:2:nargs);
    for i=1:length(names)
        if names{i} == 'cell_groupings'
            cell_groupings = values{i};
        end 
        if names{i} == 'custom_colors'
            custom_colors = values{i};
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
    write_color_tracks(custom_colors,project_directory+'color_data_gene_sets.csv')
    
    % Save gene colortracks
    disp('Writing gene colors');
    mkdir(project_directory+'gene_colors');
    II = round(length(gene_list) / 50) + 1;
    for j=0:49
        fname = project_directory+'gene_colors/color_data_all_genes-'+num2str(j)+'.csv';
        current_row = 1;
        for i = II*j:II*(j+1)
            if mean(E(:,i+II*j)) > 0.1
                all_gene_colors{current_row,1} = gene_list{current_row};
                all_gene_colors{current_row,2:} = E(:,i+II*j);
            end
        end
        write_color_tracks(all_gene_colors, fname);
    end
    
    % Save color tracks summary statistics
    disp('Saving color track statistics');
    colorstats = cell(size(E,2)+size(custom_colors,1),5)
    fieldlist = cell(size(E,2)+size(custom_colors,1))
    for i = 1:size(E,2)
        colorstats{i,1} = mean(E(:,i));
        colorstats{i,2} = std(E(:,i));
        colorstats{i,3} = 0
        colorstats{i,4} = max(E(:,i));
        colorstats{i,5} = prctile(E(:,i));
        fieldlist{i} = gene_list{i};
    end
    for i = 1:size(custom_colors,1)
        colorstats{i+size(E,2),1} = mean(custom_colors(i,:));
        colorstats{i+size(E,2),2} = std(custom_colors(i,:));
        colorstats{i+size(E,2),3} = 0
        colorstats{i+size(E,2),4} = max(custom_colors(i,:));
        colorstats{i+size(E,2),5} = prctile(custom_colors(i,:));
        fieldlist{i+size(E,2)} = custom_colors{i,1};
    end
    color_stats = cell2struct(colorstats,fieldlist,1);  
    json.write(color_stats, strcat(project_directory,'color_stats.json'));
    
	% Save categorical coloring data
    disp('Saving categorical coloring data');
    categorical_coloring_data = cell(size(cell_groupings,1));
    for i=1:size(cell_groupings,1)
        labelset = unique(cell_groupings{i,2:end});
        label_colors = cell2struct(labelset,get_categorical_colors(length(labelset)));
        categorical_coloring_data{i} = struct('label_colors',label_colors,'label_list',cell_groupings{i,2:end});
    end
    categorical_coloring_data = cell2struct(categorical_coloring_data,cell_groupings{:,1},1);
    json.write(categorical_coloring_data,strcat(project_directory+'categorical_coloring_data.json'));_
end