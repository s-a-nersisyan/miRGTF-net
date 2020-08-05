# miRGTF-net
## About
miRGTF-net is a tool for construction of miRNA-gene-TF regulatory networks using both database-level and transcriptomics data.

## Requirements
- Python 3.6+
- scikit-learn (0.23+), networkx (2.3+), pandas, numpy, scipy

## Usage
python3 run.py /path/to/config.json

All running options should be specified in config.json file (see example with default values in the repository):
- interaction_databases_path: path to folder where interaction databases are stored. A database should be presented as a two-column tab-separated (tsv) file, see examples in input_examples/interaction_databases folder (all databases listed in the manuscript)
- gene_expression_table_path: path to tsv table containing gene expression estimates. Rows of a table should correspond to genes, columns to samples, see example in input_examples/gene_expression.tsv
- miRNA_expression_table_path: same as for gene expression but for miRNAs. Note that order of samples should match for gene and miRNA expression tables (see input_examples/miRNA_expression.tsv)
- sign_constraints: any user-defined constraints on interaction directions (signs) should be put there. Do not write anything for databases without constraints. In example we assume negative correlation between miRNA and its target gene (miRTarBase_v7_0) and positive correlation between host gene and corresponding intronic miRNA
- output_path: path were all output files should be stored
- Spearman_correlation_cutoff_percentile: percentile, specifying fraction of edges which will be removed during pre-processing. For example, value equal to 90 will preserve edges with top-10% Spearman correlations. The lower value is, the more "weak" edges will appear in the constructed network
- incoming_score_threshold: threshold on regression goodness of fit (R^2 values). The higher value is, the more confident regulations will be included in the network
- interaction_score_cutoff_percentile: percentile specifying fraction of edges which should be discarded at the final step. For example, value equal to 10 will discard edges with lowest-10% interaction scores. The higher value is, the more edges will be discarded
