# here are the commands to tun the haploblock pipeline from:
# https://github.com/collaborativebioinformatics/Haploblock_Clusters_ElixirBH25

git clone https://github.com/collaborativebioinformatics/Haploblock_Clusters_ElixirBH25.git
cd Haploblock_Clusters_ElixirBH25/

# download data for chr19 as described in the README

# before running the pipeline, modify haploblock_pipeline/config/default.yaml
python3 haploblock_pipeline/main.py --config haploblock_pipeline/config/default.yaml --step 1
python3 haploblock_pipeline/main.py --config haploblock_pipeline/config/default.yaml --step 2
python3 haploblock_pipeline/main.py --config haploblock_pipeline/config/default.yaml --step 3
python3 haploblock_pipeline/main.py --config haploblock_pipeline/config/default.yaml --step 4
python3 haploblock_pipeline/main.py --config haploblock_pipeline/config/default.yaml --step 5

# analyze clusters
cd out_dir/APOE/clusters
awk '{counts[$1]++} END {for (rep in counts) print rep, counts[rep]}' chr19_44885589-44935297_cluster.tsv | sort -k2,2nr > cluster_size.tsv
less cluster_size.tsv
cat cluster_size.tsv | wc -l  # get the number of clusters
grep " 1" cluster_size.tsv | wc -l  # get the number of clusters with one individual only
