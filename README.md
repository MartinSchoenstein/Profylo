# Profylo

Profylo: an accessible phylogenetic profiling analysis python package : similarity metrics, clustering methods and modules exploration and visualisation

**Documentation: https://martinschoenstein.github.io/Profylo/**

---

## Installation

All dependencies are described in the pyproject.toml file. As well as the appropriate Python version (3.12.10 recommended) for the conda environment that will host the package.
You can create it with the following commands:
```bash
conda create -n Profylo_env python=3.12
conda activate Profylo_env
```

For now, the package is accessible by downloading the folder on GitHub. Once this is done, in the directory containing the .toml file, you can install the package with the following commands:
```bash
git clone https://github.com/MartinSchoenstein/Profylo.git
cd Profylo
pip install .
```

---

## How to use

Here are some examples of how to use Profylo, using some important functions. The data used and generated files are available on GitHub in the example folder.

Profile processing and modification functions, some used automatically by the library (such as transformation into transition vector) are available but not described here:
```python
from profylo import pre_processing as pre
```

**All functions available in the library and their different options are described in the documentation: https://martinschoenstein.github.io/Profylo/**

#### 1) Generating similarity scores (cotransition scores here) between all profiles in a matrix:

```python
from profylo import Profylo as pro

pro.distance_profiles(x = "example/profiles.csv",  tree = "example/tree.nwk", method = "cotransition", consecutive = False, path = "example/cotransition_similarity.csv")
```
Other similarity or distance metrics are available, some with their own customization options, and some without requiring a phylogenetic tree.
  

#### 2) Obtaining functional modules (using the connected components technique here):

```python
from profylo import post_processing as post 

post.graph_modules("example/cotransition_similarity.csv", distance = "cotransition", threshold = 0.3, path = "example/connected_components.txt")
```
Other clustering methods, some graph-based and non-graph-based, are available.
  

#### 3) Visualiation of a module:

- Heatmaps (It would have been possible to highlight certain clades with the *clades* argument of the function):
```python
from profylo import post_processing as post 

post.profils_heatmap("example/profiles.csv", ['P22102', 'P31939', 'Q06203', 'O15067'], tree = "example/tree.nwk", path = "example/heatmaps_cluster5.png")
```

- Annotated tree (to be viewed in external software):
```python
from profylo import post_processing as post 

post.tree_annotation(['P22102', 'P31939', 'Q06203', 'O15067'], "example/profiles.csv", path_tree = "example/tree.nwk", path = "example/annotated_tree_cluster5.nhx")
```
  

#### 4) Functional characterization of modules:

```python
from profylo import post_processing as post 

post.go_enrichment("example/connected_components.txt", gaf = "example/gaf.gaf", path = "GO_enrichment.csv")
```
  

#### 5) Obtaining phylogenetic information relating to modules, the parsimony score is used as a score of atypicality of the modules:

```python
from profylo import post_processing as post 

post.phylogenetic_statistics("example/connected_components.txt", profils = "example/profiles.csv", path_tree = "example/tree.nwk", path = "example/phylogenetic_statistics.csv")
```

---

## Crédits
Martin Schoenstein  
Pauline Mermillod  
Yannis Nevers  

CSTB - Icube 

April 2025




