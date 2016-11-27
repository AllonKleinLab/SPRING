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

SPRING is a collection of pre-processing scripts and a web browser-based tool for visualizing and interacting with high dimensional data. SPRING was developed for single cell RNA-Seq data but can be applied more generally. The minimal input is a matrix of high dimensional data points (cells) and a list of dimension names (genes). Casual users are encouraged to access our user-friendly <a href="https://kleintools.hms.harvard.edu/tools/spring.html">webserver</a>. View an example dataset <a href=https://kleintools.hms.harvard.edu/tools/springViewer.html?cgi-bin/client_datasets/centroids>here</a>. A full pyton example showing shwoing how to process your own data and boot up a local server is provided in the [Quick Start](#Quick_Start2) section. 

Low-dimensional visualizations of high-dimensional data are usually imperfect. Rather than attempting to present a single ‘definitive’ view of single cell data, SPRING allows exploration of multiple visualizations in order to develop an intuition for data structure. The core of SPRING is to create a k-nearest neighbor (kNN) graph of data points visualize the graph in 2D using a force-directed layout. A web-based interface provides a set of interactive tools to: manipulate (and thus explore) graph layout in real time; to represent any characteristic (e.g. gene expression) as a color map over the graph nodes; and to identify enriched characteristics (genes, terms) on selected graph nodes. Several export options are available to download graph representation and enriched term lists.

The SPRING subroutines can be divided into (a) pre-processing scripts that take raw inputs and convert them into data structures ready for visualization; and (b) visualization subroutines that display the pre-processed data through a web browser. The output of the pre-processing scripts is a project directory, containing a set of files with stereotyped names and formats. The visualization subroutines, implemented in javascript, accept a project directory containing the pre-processed files. We provide [pre-processing scripts](#Preprocessing2) in both Python and MATLAB. For users wishing to develop their own pre-processing scripts, a [detailed specification](#File_structures2) of the output file formats is described below. 




<a name="Installation"/>
## Installation ##

1. Download the SPRING repo: go to the green "Clone or download" button on this page
2. Alternatively: Install git and in the terminal by entering `git clone https://git@github.com/AllonKleinLab/SPRING.git`
3. If following option (2) above, you may need to change permissions using `sudo chmod -R a+w SPRING`

<a name="Quick_Start1"/>
<a name="Quick_Start2"/>
## Quick Start ##

#### Explore a pre-processed dataset ####

1. Go into the SPRING directory by entering `cd SPRING`
2. Start a local server by entering `python -m SimpleHTTPServer 8000 &`
3. In a web browser (preferably Chrome) go to <a href="http://localhost:8000/springViewer.html?datasets/centroids">http://localhost:8000/springViewer.html?datasets/centroids</a>.

#### Process your own dataset ####

_To load your own data into SPRING, they must saved to a project directory with stereotyped names and formats. We provide [preprocesing scripts](#Preprocessing3) in python and MATLAB that construct the project directory from easy inputs such as an expression matrix and a distance matrix. Sample code below uses the python preprocessing scripts to construct a project directory `datasets/frog/` from python datastructures._ 

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
## Pre-processing your data ##

SPRING pre-processing scripts and input file structures
Although simple in concept, generating a kNN graph of single cell (or other) data involves several choices. The provided subroutines provide several parameters and options for normalization, gene filtering, mesoscale dimensionality reduction, and graph construction. The routines also allow generating composite gene scores for plotting.

Input data files:

Scripts and parameters:


SPRING visualization
Instructions for launching the web interface given output of the pre-processing routines


<a name="Visualizing"/>
## Visualizing your data ##

1. Open a terminal to the SPRING directory. Create a project directory for your data if you have not done so. 
2. Start a local server by entering  `python -m SimpleHTTPServer 8000 &`
3. Go to the following URL, which must be modified with the name of your project <a href="">http://localhost:8000/springViewer.html?PATH_TO_YOUR_PROJECT_DIRECTORY</a>.
       



<a name="File_structures1"/>
<a name="File_structures2"/>
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
`GENE_NAME,cell1_expression,cell2_espression...`. For example, `Sox2,0.3,0.54,0.6... `. So if the dataset has `n` cells, this file should contain `n+1` columns. NOTE: Make sure that the file has no header, and that there are no extra commas on a line. 

2. **graph_data.json [REQUIRED]** <br>
Json file containing the graph data, with the following form (use base-0 numbering; any json compatibl format is OK):
 
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
Blah

5. **color_data_gene_sets.csv [OPTIONAL]** <br>' 
Blah

