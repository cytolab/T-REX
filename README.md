# T-REX (Tracking Responders Expanding)

### T-REX pre-print available via bioRxiv:
#### https://www.biorxiv.org/content/10.1101/2020.07.31.190454v1.full

T-REX is an algorithm developed in R for machine learning.  The T-REX acronym stands for Tracking Responders Expanding.  

T-REX takes as input two flow cytometry files that represent a pair of sampling times from one individual person.  As output, T-REX creates a UMAP analysis of an equal sampling of cells from the two files and then identifies hotspots of cells in phenotypic regions that are the most different between the two files.  

In the context of immune monitoring, the focus of T-REX on regions of great change between two sampling times can reveal disease-specific immune cells reacting to a virus.  T-REX has been tested on mass cytometry and spectral flow cytometry data from individuals with COVID-19, rhinovirus, melanoma, and leukemia.  For the associated scientific manuscript, please see Barone and Paul et al., bioRxiv 2020:

Unsupervised machine learning reveals key immune cell subsets in COVID-19, rhinovirus infection, and cancer therapy
https://www.biorxiv.org/content/10.1101/2020.07.31.190454v1.full
Sierra M. Barone*,  Alberta G.A. Paul*,  M. Muehling,  A. Lannigan,  William W. Kwok, Ronald B. Turner,  Judith A. Woodfolk,  Jonathan M. Irish

### Figure 1:

![alt text](https://www.biorxiv.org/content/biorxiv/early/2020/08/01/2020.07.31.190454/F1.large.jpg)

#### Tracking Responders EXpanding (T-REX) algorithm identifies rare cells based on significant expansion or contraction during infection or treatment.
Graphic of the Tracking Responders Expanding (T-REX) workflow. Data from paired samples of blood from a subject are collected over the course of infection and analyzed by high dimensional, high cellularity cytometry approaches (e.g., Aurora or CyTOF instrument, as with datasets here). Cells from the sample pair are then equally subsampled for UMAP analysis. A KNN search is then performed within the UMAP manifold for every cell. For every cell, the percent change between the sample pairs is calculated for the cells within its KNN region. Regions of marked expansion or contraction during infection are then analyzed to identify cell types and key features using MEM. For some datasets, additional information not used in the analysis could be assessed to determine whether identified cells were virus-specific. Finally, the average direction and magnitude of change for cells in the sample was calculated as an overall summary of how the analyzed cells changed between samples.

If you’re interested in learning more, check out the other tools on the CytoLab Github page at:
https://github.com/cytolab/

T-REX was developed for human immune monitoring using single cell cytometry in a collaboration between the laboratories of Dr. Jonathan Irish at Vanderbilt University and of Dr. Judith Woodfolk at University of Virginia.  The research was supported by the following funding resources: NIH/NCI grants U01 AI125056 (S.M.B., A.G.A.P, L.M.M., J.A.W., and J.M.I.), R01 CA226833 (J.M.I., S.M.B.), U54 CA217450 (J.M.I.), T32 AI007496 (L.M.M.) and the Vanderbilt-Ingram Cancer Center (VICC, P30 CA68485).
