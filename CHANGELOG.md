### 1.4.10

Fix cluster merging logic

### 1.4.9

Pipeline will now longer crash when all samples fail QC

### 1.4.8

Skip unchanged clusters when generating cluster sheets

### 1.4.7

Fix cluster sheet generation for new clusters

### 1.4.2

Fix merged_cluster POST

### 1.4.1

fix typo in varname

### 1.4.0

Support automated cluster merging

### 1.3.6

Fix authentication

### 1.3.4

Add API Authentication

### 1.3.3

Fix fallback tree generation

### 1.3.2

Add branch length to null tree

### 1.3.1

Fix trees

### 1.3.0

Report version to API

### 1.2.0

Deployement under CONDA_PREFIX
Adapted to API

### 0.3.2

Fixing Version
Fixing CI

### 0.3.1

Fix for the case where no cluster exists

### 0.3.0

Added distance matrices to cluster descriptions:

- distance matrix is now included in the cluster descriptions, albeit only for the main cluster as the infromation would be redundant for sub-clusters. This information is encoded as JSON records in a new field `distance_matrix`.
- the "orphan cluster" now also includes a distance matrix, including all orphans and clusters. Like this it is possible to retrieve the distance between different clusters and between orphans and clusters. THe distance to clusters uses the distance to the cluster's representative.

### 0.1.0

Initial release

