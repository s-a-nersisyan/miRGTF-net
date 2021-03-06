# miRGTF-net
## About
miRGTF-net is a tool for construction of miRNA-gene-TF regulatory networks using both database-level and transcriptomics/miRNomics data.

## Requirements
- Python 3.6+
- scikit-learn (0.23+), networkx (2.3+), pandas, numpy, scipy

## Usage
python3 run.py /path/to/config.json

All running options should be specified in config.json file (see example with default values in the repository):
- interaction_databases_path: path to folder where interaction databases are stored. A database should be presented as a two-column tab-separated (tsv) file, see examples in input_examples/interaction_databases folder (all human databases listed in the manuscript + databases for mouse)
- sign_constraints: any user-defined constraints on interaction directions (signs) should be put there. Do not write anything for databases without constraints. In example we assume negative correlation between miRNA and its target gene (miRTarBase_v8_0) and positive correlation between host gene and corresponding intronic miRNA
- gene_expression_table_path: path to tsv table containing gene expression estimates. Rows of a table should correspond to genes, columns to samples, see example in input_examples/gene_expression.tsv (100 breast cancer samples from TCGA-BRCA).
- miRNA_expression_table_path: same as for gene expression but for miRNAs. Note that order of samples should match for gene and miRNA expression tables (see input_examples/miRNA_expression.tsv)
- output_path: path were all output files should be stored
- Spearman_correlation_cutoff_percentile: percentile, specifying fraction of edges which will be removed during pre-processing. For example, value equal to 90 will preserve edges with top-10% Spearman correlations. The lower value is, the more "weak" edges will appear in the constructed network
- incoming_score_threshold: threshold on regression goodness of fit (R^2 values). The higher value is, the more confident regulations will be included in the network
- interaction_score_cutoff_percentile: percentile specifying fraction of edges which should be discarded at the final step. For example, value equal to 10 will discard edges with lowest-10% interaction scores. The higher value is, the more edges will be discarded

## TCGA data
ER+ breast cancer expression data from TCGA used in the manuscript are avaiable at https://drive.google.com/drive/folders/1c50BIlfF2SlnswJ2iTD0pXG6Rzt0Rs0m?usp=sharing.

For user convenience, we also added instruction and prepared a script how to import arbitrary TCGA data subset:
1. Access GDC Data Portal (https://portal.gdc.cancer.gov/)
2. Select a project (e.g. TCGA-BRCA), add *.FPKM.txt.gz (RNA-seq) and *.isoforms.quantification.txt (miRNA-seq) files to the cart
3. Download the cart and Sample Sheet, unpack *.tar.gz archive
4. Go to the TCGA folder and run python3 process.py /path/to/sample_sheet.tsv /path/to/gdc_download_unpacked_folder . This will generate gene and miRNA expression tables (TPM) for matched primary tumors as described in the manuscript. GENCODE v22 and miRBase v22.1 are used as default versions.

## Note about miRNA/gene identifiers
Please be careful when using different miRNA/gene/TF databases – they could have built on different versions of miRBase/RefSeq/Ensembl, which will result in missing nodes/edges. We recommend using miRBaseConverter (https://taoshengxu.shinyapps.io/mirbaseconverter/) and HGNC Multi-symbol checker (https://www.genenames.org/tools/multi-symbol-checker/) for miRNA/gene name convertations. All pre-included miRNA databases (both for human and mouse) are already compatible with miRBase v22.1.
