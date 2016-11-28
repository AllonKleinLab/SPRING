# SPRING

#### Table of Contents  
[Overview](#Overview)  
[Installation](#Installation)  
[Quick Start](#Quick_Start1)   
[Pre-processing your data](#Preprocessing1)   
[Visualizing your data](#Visualizing)   
[SPRING file structures](#File_structures1)   


<a name="Overview"/>
## Overview ##

SPRING is a collection of pre-processing scripts and a web browser-based tool for visualizing and interacting with high dimensional data. View an example dataset <a href=https://kleintools.hms.harvard.edu/tools/springViewer.html?cgi-bin/client_datasets/centroids>here</a>. SPRING was developed for single cell RNA-Seq data but can be applied more generally. The minimal input is a matrix of high dimensional data points (cells) and a list of dimension names (genes). Casual users are encouraged to access our user-friendly <a href="https://kleintools.hms.harvard.edu/tools/spring.html">webserver</a>. Heavy users and those wanting more control over the data processing pipeline may use the local installation (i.e. this github repo). A full python example showing shwoing how to process your own data and boot up a local server is provided in the [Quick Start](#Quick_Start2) section. 

Low-dimensional visualizations of high-dimensional data are usually imperfect. Rather than attempting to present a single ‘definitive’ view of single cell data, SPRING allows exploration of multiple visualizations in order to develop an intuition for data structure. The core of SPRING is to create a k-nearest neighbor (kNN) graph of data points and visualize the graph in 2D using a force-directed layout. A web-based interface provides a set of interactive tools to: manipulate (and thus explore) graph layout in real time; represent any characteristic (e.g. gene expression) as a color map over the graph nodes; and identify enriched characteristics (genes, terms) on selected graph nodes. Several export options are available to download the graph representation and enriched term lists.

The SPRING subroutines can be divided into (a) pre-processing scripts that take raw inputs and convert them into data structures ready for visualization; and (b) visualization subroutines that display the pre-processed data through a web browser. The output of the pre-processing scripts is a project directory, containing a set of files with stereotyped names and formats. The visualization subroutines, implemented in javascript, accept a project directory containing the pre-processed files. We provide [pre-processing scripts](#Preprocessing2) in both Python and MATLAB. For users wishing to develop their own pre-processing scripts, a [detailed specification](#File_structures2) of the output file formats is described below. 




<a name="Installation"/>
## Installation ##

1. Download the SPRING repo: go to the green "Clone or download" button on this page
2. Alternatively: Make sure git is installed and in the terminal enter <br> `git clone https://git@github.com/AllonKleinLab/SPRING.git`
3. If following option (2) above, you may need to change permissions using `sudo chmod -R a+w SPRING`

<a name="Quick_Start1"/>
<a name="Quick_Start2"/>
<a name="Quick_Start3"/>
## Quick Start ##

#### Explore pre-processed dataset using a local webserver ####

1. Go into the SPRING directory by entering `cd SPRING`
2. Start a local server by entering `python -m SimpleHTTPServer 8000 &`
3. In a web browser (preferably Chrome) go to <a href="http://localhost:8000/springViewer.html?datasets/centroids">http://localhost:8000/springViewer.html?datasets/centroids</a>.

#### Process your own dataset ####

_To load your own data into SPRING, the data must saved to a project directory as files with stereotyped names and formats. We provide [preprocesing scripts](#Preprocessing3) in python and MATLAB that construct the project directory from easy inputs such as an expression matrix and a distance matrix. Sample code below uses the python preprocessing scripts to construct a project directory `datasets/frog/` from python datastructures._ 

1. Unzip `example_inputs/E.npy.zip`
2. In the SPRING directory, run the following python code.

            import pickle, numpy as np

            # Import SPRING helper functions
            from preprocesing_python import *

            # Import expression matrix; rows are cells and columns are genes
            ### ****** Make sure E.npy is unzipped *************
            print 'Loading expression matrix'
            E = np.load('example_inputs/python_E.npy')

            # Filter out cells with fewer than 1000 UMIs
            print 'Filtering cells'
            E,cell_filter = filter_cells(E,1000)

            # Normalize gene expression data
            # Only use genes that make up <
            print 'Row-normalizing'
            E = row_normalize(E)

            # Filter genes with mean expression < 0.1 and fano factor < 3
            print 'Filtering genes'
            _,gene_filter = filter_genes(E,0.1,3)

            # Z-score the gene-filtered expression matrix and do PCA with 20 pcs
            print 'Zscoring and PCA'
            Epca = get_PCA(Zscore(E[:,gene_filter]),20)

            # get euclidean distances in the PC space
            print 'Getting distance matrix'
            D = get_distance_matrix(Epca)

            # load additional data
            # gene_list is a list of genes with length E.shape[1]
            # cell_groupings is a dict of the form: { <grouping_name> : [<cell1_label>, <cell2_label>,...] }
            # a "grouping" could be the sample id, cluster label, or any other categorical variable
            gene_list, cell_groupings = pickle.load(open('example_inputs/python_data.p'))

            # save a SPRING plots with k=5 edges per node in the directory "datasets/frog/"
            print 'Saving SPRING plot'
            save_spring_dir(E,D,5,gene_list,'datasets/frog', cell_groupings=cell_groupings)

3. If you haven't already, start a local server by entering `python -m SimpleHTTPServer 8000 &`
4. In a web browser, go to <a href="http://localhost:8000/springViewer.html?datasets/frog">http://localhost:8000/springViewer.html?datasets/frog</a>.


<a name="Preprocessing1"/>
<a name="Preprocessing2"/>
<a name="Preprocessing3"/>
<a name="Preprocessing4"/>
## Pre-processing your data ##

We provide pre-processing scripts in python and MATLAB that help process basic inputs into the [special files](#File_structures3) that are read by SPRING. The main function, `save_spring_dir` actually writes the project directory, taking an expression matrix and pairwaise distance matrix as inputs. The remaining functions implement basic filtering and normaization routines to produce the required distance matrix. 

### Python preprocessing ###

A full example running python pre-processing functions on example inputs is provided in the [Quick Start](#Quick_Start3) section.

### MATLAB preprocessing ###

The following code snippet will begin with basic MATLAB data structures and use them create the project directory `datasets/frog/`. To run the code, open MATLAB and go to the directory `SPRING/preprocessing_matlab/`

            % Load example impit data
            % This loads [E, gene_list, custom_colors, cell_groupings]
            % "E"              array of gene expression values (each row corresponds to a cell and each column to a gene)
            % "gene_list"      cell array of gene names with length size(E,2)
            % "custom_colors"  cell array with one row for each custom color track. The first entry in each row is the 
            %                  name of the track and subsequent entries are values for each cell. So if there are T custom 
            %                  color tracks and N cells, this should be a T x (N+1) cell array. 
            % "cell groupings" cell array with one row for each cell grouping. The first entry in each row is the name of the
            %                  of the grouping (e.g. "sampleID") and each subsequet entry is a cell label (e.g. "sample_1").
            %                  If there are T different groupungs and N cells, then this should be a T x (N+1) cell array. 
            load('../example_inputs/matlab_data.m','-mat');

            % Make sure all genes, cell groupings and custom color names can be used as fields in a struct
            % That means they cannot begin with a digit or contain "-", ".", " ", or "/"
            gene_list      = struct_field_qualified(gene_list);
            cell_groupings = struct_field_qualified(cell_groupings);
            custom_colors  = struct_field_qualified(custom_colors);

            % Filter out cells with fewer than 1000 UMIs
            disp('Filtering cells');
            [E,cell_filter] = filter_cells(E,1000);

            % Normalize gene expression data
            % Only use genes making up <5% of total UMIs
            disp('Row-normalizing');
            E = row_normalize(E);

            % Filter genes with mean expression < 0.1 and fano factor < 3
            disp('Filtering genes');
            [~,gene_filter] = filter_genes(E,0.1,3);

            % Z-score the gene-filtered expression matrix and do PCA with 20 pcs
            disp('Zscoring and PCA');
            [coeff,score,latent] = pca(zscore(E(:,gene_filter)));
            Epca = score(:,1:20); 

            % get euclidean distances in the PC space
            disp('Getting distance matrix');
            D = pdist2(Epca,Epca);

            % save a SPRING plots with k=5 edges per node in the directory "../datasets/frog/"
            disp('Saving SPRING plot');
            save_spring_dir(E,D,5,gene_list,'../datasets/frog2', 'cell_groupings',cell_groupings,'custom_colors',custom_colors);

<a name="Visualizing"/>
## Visualizing your data ##
0. At this point, it is assumed that you have already created a project directory using the pre-processing [scripts](#Preprocessing4). 
1. Open a terminal to the SPRING directory. 
2. Start a local server by entering  `python -m SimpleHTTPServer 8000 &`
3. Go to the following URL, which must be modified with the name of your project <a href="">http://localhost:8000/springViewer.html?PATH_TO_YOUR_PROJECT_DIRECTORY</a>.
       



<a name="File_structures1"/>
<a name="File_structures2"/>
<a name="File_structures3"/>
## SPRING file structures ##

The SPRING project directory must contain files with stereotyped names and formats. Matlab and Python [scripts](#Preprocessing3) are provided to create these files. Here is a guide to the file names and formats:
 
1. **gene_colors/color_data_all_genes-*.csv [REQUIRED]** <br>
In a directory called `gene_colors` there must be (at most 50) base-0 numbered files called `color_data_all_genes-*.csv.` e.g.
 
            color_data_all_genes-0.csv
            color_data_all_genes-1.csv
            color_data_all_genes-2.csv
            ...
            color_data_all_genes-50.csv
Each of these files should contain gene expression for a subset of genes, with one gene on each row. The rows have the following format:
`GENE_NAME,cell1_expression,cell2_espression...`. For example, `Sox2,0.3,0.54,0.6... `. So if the dataset has `n` cells, this file should contain `n+1` columns. NOTE: Make sure that the file has no header. 

2. **graph_data.json [REQUIRED]** <br>
Json file containing the graph data, with the following form (use base-0 numbering; any json compatible format is OK):
 
            {  "nodes": [ {   "name": cell0, "number": 0 },     // List of nodes
                          {   "name": cell1, "number": 1 },
                          ....
                          {   "name": cellN, "number": N } ],
               "links": [ { "source": 10, "target": 23 },      // List of edges
                          { "source": 29, "target": 50 },
                          ....
                          { "source": 40, "target": 125 }  ] }

3. **color_stats.json [REQUIRED]** <br>'
Json file containing pre-calculated summary statistics of the various coloring tracks, including those in `color_data_gene_sets.csv` (see below) and `color_data_all_genes-*.csv`. The file should contain a dictionary mapping each color track-name to a list  [MEAN, STANDARD DEVIATION, MIN, MAX, 99-PERCENTILE] with summary statistics for that color track. Thus, this file could have the form:

            {  "Sox2":  [ 0.1, 0.2, 0, 1.46, 1.22],
               "Brca":  [ 5.2, 4.1, 0, 20.3, 18.1],
               ...
               "Gata1": [ 0.4, 0.3, 0, 5.42, 4.11]  }

4. **categorical_coloring_data.json [OPTIONAL]** <br>' 
Json file containing cell groupings, i.e. categorical variables such as sample ID, cluster label, etc. For each cell grouping, a color map and label list must be provided, as follows.

            { "SampleID":  {  "label_list": [  "Sample1", "Sample2", "Sample2", ... "Sample1"  ], 
                              "label_colors": { "Sample1": "#00007f", 
                                                "Sample2": "#00007f" }  },

              "ClusterID": {  "label_list": [  "Cluster1", "Cluster3", ... "Cluster2"  ], 
                              "label_colors": { "Cluster1": "#00007f", 
                                                "Cluster2": "#00007f",
                                                "Cluster3": "#00007f"}  }  }                                   
The "label_list" array should contain one string for each cell. The "label_colors" map should have one name-color pair for each distinct cell label in "label_list".                               


5. **color_data_gene_sets.csv [OPTIONAL]** <br>
This csv file stores continuous variables for coloring the data, such as signature scores or cell pseudotime. Each line of the file corresponds to one coloring track, whith the name of the track followed by a sequence of values: `TRACK_NAME,cell1_value,cell2_value...`. For example, `Cell_cycle_score,0.3,0.54,0.2... ` So if the dataset has `n` cells, this file should contain `n+1` columns. NOTE: Make sure that the file has no header. 

