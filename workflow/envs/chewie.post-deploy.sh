#!/usr/bin/env bash
set -Eeu

# Repo URL
chewie_repo="https://gitlab.com/bfr_bioinformatics/chewieSnake.git"

# Commit hash to use
commit="b7117542d3460c8234ab5a3a4f5545aef443cbf2"

# Local directory to save the Repo
local_dir="$HOME/.nrw-geuebt/chewieSnake/"

# if already exists, wipe it clean and redo clone
[ -d "$local_dir" ] && rm -rf "$local_dir"

# clone and checkout repo
echo "Cloning ChewieSnake and checking out stable commit"
git clone -q "$chewie_repo" "$local_dir"
cd "$local_dir"
git checkout "$commit"

# Patch wrapper and Snakefile
echo "Applying patches"
patch -s --directory="$local_dir" --strip=1 << END
diff --unified --recursive --no-dereference chewieSnake-orig/chewieSnake_join.py chewieSnake/chewieSnake_join.py
--- chewieSnake-orig/chewieSnake_join.py	2023-05-30 10:12:16.855136923 +0200
+++ chewieSnake/chewieSnake_join.py	2023-05-30 10:14:42.000000000 +0200
@@ -147,7 +147,7 @@
     parser.add_argument('--clustering_method', help = 'The agglomeration method to be used for hierarchical clustering. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC); default: single',
                         default = "single", required = False)
     parser.add_argument('--distance_threshold', help = 'A single distance threshold for the exctraction of sub-trees; default: 10',
-                        default = 10, required = False)
+                        default = 10, type=int, required = False)
     parser.add_argument('--cluster_names', help = 'A file with potential names for cluster names, one name per line, default: repo_path/scripts/cluster_names_reservoir.txt',
                         default = os.path.join(repo_path, "scripts", "cluster_names_reservoir.txt"), required = False)
     parser.add_argument('--subcluster_thresholds', help = 'A list of distance thresholds for subclustering; default: [3]',
diff --unified --recursive --no-dereference chewieSnake-orig/chewieSnake_join.smk chewieSnake/chewieSnake_join.smk
--- chewieSnake-orig/chewieSnake_join.smk	2023-05-30 10:12:16.451155752 +0200
+++ chewieSnake/chewieSnake_join.smk	2023-05-30 10:14:52.000000000 +0200
@@ -76,7 +76,7 @@
         "merged_db/sample_cluster_information.tsv", # from rule clustering
         #"merged_db/samples2clusters.csv" # result from checkpoint
         #"merged_db/listofreports.txt", # result from checkpoint
-        "merged_db/report.html"
+        # "merged_db/report.html"
 
 rule all_reportonly:
     input:
diff --unified --recursive --no-dereference chewieSnake-orig/envs/grapetree3.yaml chewieSnake/envs/grapetree3.yaml
--- chewieSnake-orig/envs/grapetree3.yaml	2023-05-30 10:17:00.031148266 +0200
+++ chewieSnake/envs/grapetree3.yaml	2023-06-01 12:08:18.000000000 +0200
@@ -3,50 +3,5 @@
   - conda-forge
   - bioconda
   - defaults
-  - r
 dependencies:
-  - _libgcc_mutex=0.1
-  - _openmp_mutex=4.5
-  - ca-certificates=2020.11.8
-  - certifi=2020.11.8
-  - ld_impl_linux-64=2.35.1
-  - libffi=3.2.1
-  - libgcc-ng=9.3.0
-  - libgomp=9.3.0
-  - libstdcxx-ng=9.3.0
-  - ncurses=6.2
-  - openssl=1.1.1h
-  - pip=20.2.4
-  - python=3.8.6
-  - python_abi=3.8
-  - readline=8.0
-  - sqlite=3.33.0
-  - tk=8.6.10
-  - wheel=0.35.1
-  - xz=5.2.5
-  - zlib=1.2.11
-  - pip:
-    - chardet==3.0.4
-    - click==7.1.2
-    - decorator==4.4.2
-    - ete3==3.1.2
-    - flask==1.1.2
-    - grapetree==2.2
-    - idna==2.10
-    - itsdangerous==1.1.0
-    - jinja2==2.11.2
-    - llvmlite==0.34.0
-    - markupsafe==1.1.1
-    - networkx==2.5
-    - numba==0.51.2
-    - numpy==1.19.4
-    - pandas==1.1.4
-    - psutil==5.7.3
-    - python-dateutil==2.8.1
-    - pytz==2020.4
-    - requests==2.25.0
-    - setuptools==50.3.2
-    - six==1.15.0
-    - unidecode==1.1.1
-    - urllib3==1.26.2
-    - werkzeug==1.0.1
+  - grapetree=2.1

END

