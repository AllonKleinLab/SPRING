# SPRING

SPRING is a tool for exploring topological relationships in single-cell sequencing data that relies on interactive visualization of a force-directed graph. In the SPRING graph, each node is a cell and each edge is an instance of proximity in gene expression space. Users ca upload their own expression data to SPRING using our webserver https://kleintools.hms.harvard.edu/tools/spring.html

However, heavy users are encouraged to create a local installation of SPRING by following these steps:

0. Install git
1. Clone the SPRING repo from github. In the terminal enter `git clone git@github.com:AllonKleinLab/SPRING.git`

To run SPRING locally

2. Go to the SPRING directory by entering `cd SPRING`
3. Start a local server by entering `python -m SimpleHTTPServer 8000 &`
4. Open a browser (Chrome is recommended) and go to http://localhost:8000/springViewer.html?datasets/centroids

To learn more, see <a href="http://en.wikipedia.org/wiki/Main_Page">Wikipedia</a>.
