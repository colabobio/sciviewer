# Single Cell Interactive Viewer (Sciviewer)

This is an interactive viewer for 2D embeddings such as [UMAP](https://umap-learn.readthedocs.io/en/latest/) or [tSNE](https://lvdmaaten.github.io/tsne/) of high dimensional single-cell RNA-seq data that is run directly out of the [Jupyter Notebook](https://jupyter.org/) environment. The user can select cells in the 2D plane and the viewer will calculate the differential expression between the selected and the unselected cells. Alternatively, the user can select a group of cells and a direction and the viewer will identify the genes with the greatest variation (Pearson correlation) along that direction. See a video of how this works below. Also, see the example tutorial in this repository (a small example of [3000 PBMCs](tutorials/sciviewer_example_3K_PBMC.ipynb) that illustrates all the input options, and a larger one of [50,000 circulating T-cells](tutorials/sciviewer_example_50K_Tcell.ipynb)).

***Important note:*** The code here is an initial proof-of-concept, the development of Sciviewer continues by extending [anywidget](https://anywidget.dev/), [regl-scatterplot](https://github.com/colabobio/regl-scatterplot), and [jupyter-scatter](https://github.com/colabobio/jupyter-scatter) to recreate it as a fully-featured visualization tool for Jupyter notebooks and completely integrated into the [scverse](https://scverse.org/) ecosystem:

![image](https://github.com/user-attachments/assets/59c2da50-7e09-4a14-a42c-2facba8c1779)
*Mockup of the next version of Sciviewer based on anywidget, regl-scatterplot, and jupyter-scatter*

## Installation

The main requirement for sciviewer is [py5](http://py5.ixora.io/install/) which in turn requires Python 3.8. We recommend using the conda package manager to install the necessary dependencies fo sciviewer. Conda can be installed following the instructions [here](https://docs.conda.io/en/latest/miniconda.html). Then follow the steps below to install sciviewer.

1. Prepare and activate the conda environment containing dependencies for py5:

```
conda env create -n sciviewer -f https://raw.githubusercontent.com/colabobio/sciviewer/master/sciviewer-env.yml
conda activate sciviewer
```
Alternatively, if you want to append the needed dependencies to an existing conda environment, instead of creating a new one, you can do the following:

```
conda env update -n your_existing_environment -f https://raw.githubusercontent.com/colabobio/sciviewer/master/sciviewer-env.yml
conda activate your_existing_environment
```

Note, Sciviewer currently requires Python 3.8 or greater.

2. Install OpenJDK 1.7. Py5 does not work with OpenJDK 11, which is the one available through Conda at the moment. After creating an activating the sciviewer environment, there are two options:

* If you already have OpenJDK 17 installed in your system, you can make it available to sciviewer by setting the JAVA_HOME environmental variable. For instance, if you installed [Adoptium OpenJDK 17](Adoptium), the home folder should be ```/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home``` so you can get the JAVA_HOME as follows:

```
 export JAVA_HOME=/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home
```

* Otherwise, you can install OpenJDK 17 as explained in Py5's [Anaconda setup](http://py5.ixora.io/content/install.html#brief-steps).


3. Next, install the sciviewer package using pip:
```
pip install sciviewer
```

And that is it, the module is now installed and ready to be used.

To uninstall, use:
```
pip uninstall sciviewer
```

3. Now launch jupyter from within the activated conda environment and you are good to go.

```
jupyter lab
```

## Quick start

Sciviewer is executed from a Jupyter notebook such as in the examples directory. It is run by initializing a SCIViewer object with the 2D embedding (# cells X 2) and the expression data (# cells X # genes) and then running the explore_data method. E.g.

```
from sciviewer import SCIViewer
svobj = SCIViewer(expr, umap)
svobj.explore_data()
```

Running the code above will cause the visualizer to appear. The umap and expression data can also now be provided directly as a [Scanpy AnnData](https://scanpy.readthedocs.io/en/stable/usage-principles.html#anndata), see below or [tutorials](./tutorials/sciviewer_example_3K_PBMC.ipynb) for details.

Click the video link below for a ~3 minute tutorial on how to use the visualizer:
[![Watch the video](https://img.youtube.com/vi/YgvMmvgFFE0/maxresdefault.jpg)](https://youtu.be/YgvMmvgFFE0)

Note, if you are running Sciviewer from the Jupyter notebook as in the tutorials, you need to load the py5 magic extension before you can start the visualizer

```
%load_ext py5
```

and if you are running it from a jupyter notebook on a mac computer, you need to add an additional magic extension

```
%gui osx
```

## Some key usage points:
 - __Inputs__: The expression data can be provided as a Scanpy AnnData, Pandas DataFrame, a Numpy ndarray, or as a scipy sparse [csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html).
 - __AnnData expression__: For AnnData objects, the expression data are accessed from the .X attribute by default. Setting the use_raw argument to True causes it to be accessed from the .raw.X attribute instead. If the data are sparse, sciviewer requires it to be in the csc_matrix format. See the tutorial for how to convert between sparse matrix formats
 - __Sparsity__: Providing the data as a sparse csc_matrix is recommended for large datasets as it can lead to a considerable (1-2 order or magnitude) performance speedup. See [this notebook](./tutorials/sciviewer_example_3K_PBMC.ipynb) as an example.
 - __Gene/cell names__: If the expression data is provided as a Pandas DataFrame, the cell names are inferred from the index and the gene names are inferred from the columns. If it is provided as a Scanpy AnnData, the gene names come from the index of the .var attribute and the cell names come from the index of the .obs attribute. Otherwise, the gene names and cell names can be provided when initializing the SCIViewer class with the gene_names and cell_names arguments, or will be initialized with generic names. 
 - __Real time updating of python variables__ The selected_cells attribute of the sciviewer object is updated whenever a new set of cells are selected, regardless of the mode, and contains information about the selected cells. The results_proj_correlation attribute of the sciviewer object is updated whenever a new selection is made in the "directional" mode and contains the Pearson correlation and P-values of all genes for the selected direction and cells. The results_diffexpr attribute is updated when a new selection is made in the "differential" mode and contains the T-statistic and P-value for the differential expression test (simple Welch's T-test). These are updated in real time as the visualizer is in use.

See the tutorial notebooks for more details

## Development / debugging

For development purposes, it can be helpful to import sciviewer directly rather than installing the package. See the extras/debugging directory for notebooks with examples of how to do this e.g. [debug_example_3K_PBMC.ipynb](extras/debugging/debug_example_3K_PBMC.ipynb).

## Publications

* Kotliar D, Colubri A. Sciviewer enables interactive visual interrogation of single-cell RNA-Seq data from the python programming environment. (2021). Bioinformatics. doi: [10.1093/bioinformatics/btab689](https://academic.oup.com/bioinformatics/article/37/21/3961/6380550)
* Kotliar D, Colubri A. Sciviewer enables interactive visual interrogation of single-cell RNA-Seq data from the Python programming environment. (2021). bioRxiv. doi: [10.1101/2021.08.12.455997](https://www.biorxiv.org/content/10.1101/2021.08.12.455997v1) (preprint)
