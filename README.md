# Single Cell Interactive Viewer

This is an interactive viewer for Jupyter notebooks that allows users exploring an embedding of single-cell RNA-seq data. The user can select directions in the 2D embedding plane (e.g.: UMAP or tSNE), and the viewer will calculate the genes with the highest correlation along those directions. In addition, you can select a group of cells and find the genes that are most differentially expressed.

## Installation

- Prepare conda environment containing dependencies for py5:
  ```conda env create -n sciviewer -f https://raw.githubusercontent.com/colabobio/sciviewer/master/sciviewer-env.yml```
- Activate the environment created above:
  ```conda activate sciviewer```
- Alternatively, if you want to append the needed dependencies to an existing conda environment, you can do the following:
   -  ```conda env update -n your_existing_environment -f https://raw.githubusercontent.com/colabobio/sciviewer/master/sciviewer-env.yml```
   -  ```conda activate your_existing_environment```
- Install the sciviewer package:
  ```python -m pip install --index-url https://test.pypi.org/simple/ --no-deps sciviewer```
- To uninstall:
  ```pip uninstall sciviewer```
- To load the module from a notebook wihtout installing (good for debugging), add the following imports to the notebook:

```
 import sys
 sys.path.insert(0, '../sciviewer/')
 from sciviewer import SCIViewer
  ```
- Note, the above won't work if you don't remove the relative imports from the local .py files.

## Jupyter notebook

This Jupyter notebook demonstrates the use of [py5](http://py5.ixora.io/), a version of Processing for Python, to create an interactive viewer of RNA-seq embedding data. Find it under the notebook folder.

Processing's drawing API and engine can be run from a Jupyter notebook to enable users explore their datasets and then bring the results of their exploration directly into the notebook for further analysis:

![RNA-seq embed viewer running a Jupyter notebook](extras/images/jupyter_screenshot1.jpg)

The datasets can be quite large, as seen in the next image:

![RNA-seq embed viewer loading a large dataset](extras/images/jupyter_screenshot2.jpg)
