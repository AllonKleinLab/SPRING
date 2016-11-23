# SPRING

SPRING is a tool for exploring topological relationships in single-cell sequencing data that relies on interactive visualization of a force-directed graph. In the SPRING graph, each node is a cell and each edge is an instance of proximity in gene expression space. Users ca upload their own expression data to SPRING using our webserver https://kleintools.hms.harvard.edu/tools/spring.html

However, heavy users are encouraged to create a local installation of SPRING by following these steps:

1. Install git
2. Clone the SPRING repo from github. In the terminal enter `git clone git@github.com:AllonKleinLab/SPRING.git`

SPRING is now installed. Next, load a test dataset by:

3. Go to the SPRING directory by entering `cd SPRING`
4. Start a local server by entering `python -m SimpleHTTPServer 8000 &`
5. In a web browser (preferably Chrome) go to <a href="http://localhost:8000/springViewer.html?datasets/centroids">http://localhost:8000/springViewer.html?datasets/centroids</a>.
