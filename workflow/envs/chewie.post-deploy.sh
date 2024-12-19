#!/usr/bin/env bash
set -Eeu

# Repo URL
chewie_repo="https://gitlab.com/bfr_bioinformatics/chewieSnake.git"

# Commit hash to use
commit="b7117542d3460c8234ab5a3a4f5545aef443cbf2"

VERSION=$(cat "../../VERSION")

# Local directory to save the Repo
local_dir="$HOME/.nrw-geuebt/geuebt-chewie-${VERSION}/chewieSnake/"

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

END

