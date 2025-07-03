TMT Proteomics Analysis Pipeline
================================

The TMT pipeline provides comprehensive proteomics analysis with advanced normalization, batch correction, and differential expression using Limma.

Script Location
---------------

``DE_pipelines/TMT/tmt_command_line_script.py``

Input Data Requirements
-----------------------

Expression Data
~~~~~~~~~~~~~~

* **Required column**: ``ProteinID`` containing protein identifiers
* **Sample columns**: Additional columns for each sample, where rows represent protein expression levels
* **Format**: CSV file

Protein Metadata
~~~~~~~~~~~~~~~

* **Required columns**: ``ProteinID`` (same as in expression data), ``GeneSymbol`` (corresponding gene names)
* **Optional columns**: 
  
  * ``Peptides`` (number of unique peptides, required if peptide normalization requested)
  * ``Sequence`` (required for canonical isoform identification using UniParc)
  
* **Peptide format**: Information on peptides used for detection in format: ``"{'Peptide_A','Peptide_B'}"``
* **Format**: CSV file

Sample Metadata
~~~~~~~~~~~~~~

* **Required columns**: 
  
  * ``TMT ID`` specifying sample names
  * ``Condition`` specifying biological conditions
  * ``Replicate`` specifying replicates
  * Both columns are expected to contain string values

* **Additional columns**: batch information, covariates to preserve biological information during batch correction
* **Format**: Comma delimited file

Complete Parameter Reference
----------------------------

Basic Arguments
~~~~~~~~~~~~~~

``--analysis_type {preprocessing,de,both}``
    Type of analysis to run: preprocessing, differential expression (de), or both

``--organism {homo_sapiens,mus_musculus,other}``
    Organism to use for reference files and databases

``--preprocessing_config PREPROCESSING_CONFIG``
    Path to JSON preprocessing configuration file

``--de_config DE_CONFIG``
    Path to JSON de configuration file

``--verbose``
    Enable verbose output

``--ignore_excel_limit``
    Boolean flag to save all results (for differential testing or enrichment) together regardless if number of combined rows exceeds excel display limit. (default False)

Preprocessing Parameters
~~~~~~~~~~~~~~~~~~~~~~~

``--expression_data_path EXPRESSION_DATA_PATH``
    Path to the expression data file. The file must include a "ProteinID" column containing protein identifiers and additional columns for each sample, where rows represent protein expression levels.

``--protein_metadata_path PROTEIN_METADATA_PATH``
    Path to the protein metadata file. The file must include ``ProteinID`` column same as in expression data file. ``GeneSymbol`` should contain corresponding gene names for each protein. ``Peptides`` column with number of unique peptides used for detection should be present in the file, if peptide normalization is requested. If user wants to identify canonical isoforms using UniParc, ``Sequence`` column must be present. Rows must contain information on peptides used for detection of a given protein in following format. Ex. ("{{'Peptide_A','Peptide_B'}}")

``--sample_metadata_path SAMPLE_METADATA_PATH``
    Path to the comma delimited sample metadata file. The file must include 'TMT ID' column specifying sample names , 'Condition' columns specifying biological conditions and column 'Replicate' specifying replicates. 'TMT ID' and 'Replicate' columns are expected to contain string values.

``--preprocessing_output_dir PREPROCESSING_OUTPUT_DIR``
    Directory where preprocessing results will be saved

``--h5ad_filename H5AD_FILENAME``
    Filename for the h5ad output file. This file will be used for downstream differential expression analysis (default: proteomics_data.h5ad)

``--mitocarta_path MITOCARTA_PATH``
    Path to the MitoCarta data to identify genes related OXPHOS complexes (default: species-specific, for example /omics_downstream_pipeline/ref_files/human/human_genes_mitocarta3.0.csv.gz)

``--uniparc_db_path UNIPARC_DB_PATH``
    Path to UniParc database (required if ``--identify_canonical_proteins`` is True and ``--apply_uniparc_mapping`` is specified). Default file contains information on homo sapiens and mus musculus proteins (default: None)

``--swiss_prot_ref_path SWISS_PROT_REF_PATH``
    Path to Swiss-Prot reference file that would be used for canonical mapping or proteins (required if identify_canonical_proteins is True). Default file contains information on homo sapiens and mus musculus proteins (default: None)

``--min_mean MIN_MEAN``
    Minimum mean expression threshold used for identification of highly variable proteins by scanpy.pp.highly_variable_genes (default: 0.5)

``--max_mean MAX_MEAN``
    Maximum mean expression threshold used for identification of highly variable proteins by scanpy.pp.highly_variable_genes (default: 8)

``--min_disp MIN_DISP``
    Minimum dispersion threshold used for identification of highly variable proteins by scanpy.pp.highly_variable_genes (default: 1)

``--skip_sum_normalize``
    Flag to skip normalizing data by per sample sum (default: False)

``--target_sum TARGET_SUM``
    Target sum per sample, used if skip_sum_normalize is set to False (default: 1000000.0)

``--skip_peptide_normalize``
    Flag to skip normalizing by peptide count using ``Peptides`` column. Column should be present, if peptide normalization is requested.(default: False)

``--skip_mad_normalization``
    Flag to skip Median Absolute Deviation normalization of the data (default: False)

``--batch_correct``
    Flag whether to perform batch correction. If specified, column specified in batch_id_column argument will be used for batch correction using Combat package. Input layer for correction is dependent on normalization -related flags ``--skip_sum_normalize``, ``--skip_peptide_normalize``, ``--skip_mad_normalization``. Most normalized layer is used for correction. (default: False)

``--batch_id_column BATCH_ID_COLUMN``
    Batch information column in anndata_obj.obs that will be used for correction. Batches in column should be numerically encoded with integers. Used if batch_correct is True (default: None)

``--mod_id_column MOD_ID_COLUMN``
    Column name in sample metadata for outcome of interest and other covariates besides batch. Used to create model matrix to preserve biological variation during correction. Used if ``--batch_correct`` is specified.(default: None)

``--identify_canonical_proteins``
    Flag to identify canonical protein isoforms. If True, ``--swiss_prot_ref_path`` should be specified (default: False)

``--apply_uniparc_mapping``
    Use UniParc database for the most up-to-date ProteinID. ``--uniparc_db_path`` should be specified (default: False)

``--scale_max_value SCALE_MAX_VALUE``
    Maximum scaling value for scanpy.pp.scale (default: 10)

``--pca_color PCA_COLOR [PCA_COLOR ...]``
    List of column names to use as color categories in PCA plots (default: ['Condition', 'Replicate'])

``--pca_components PCA_COMPONENTS [PCA_COMPONENTS ...]``
    List of comma-separated pairs of PCA components (default: ['1,2', '2,3', '1,3'])

``--umap_plot``
    Generate UMAP plots (default: False)

``--n_neighbors N_NEIGHBORS``
    Number of nearest neighbors to compute UMAP on (default: 5)

``--n_pcs N_PCS``
    Number of principal components to use for UMAP construction (default: 9)

``--skip_plot_oxphos``
    Flag to skip plotting OXPHOS heatmaps (default: False)

``--oxphos_complexes OXPHOS_COMPLEXES [OXPHOS_COMPLEXES ...]``
    List of strings of OXPHOS complexes names to generate heatmaps for (default: ['Complex I', 'Complex II', 'Complex III', 'Complex IV', 'Complex V'])

Differential Expression Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--input_file INPUT_FILE``
    Path to the AnnData file for DE analysis

``--de_output_dir DE_OUTPUT_DIR``
    Directory where DE results will be saved

``--anndata_data_layer ANNDATA_DATA_LAYER``
    AnnData layer to use for DE analysis (default: mad_log1p_sum_norm_peptide_norm)

``--ignore_non_canonical``
    Flag to ignore non canonical values in DE analysis (default: False)

``--design_factors DESIGN_FACTORS [DESIGN_FACTORS ...]``
    Sample metadata columns to include as design factors for Limma differential analysis.(default ['Condition'])

``--condition_pairs CONDITION_PAIRS [CONDITION_PAIRS ...]``
    Condition pairs for comparison (format: cond1 cond2 [cond3 cond4 ...])

``--logfc_threshold LOGFC_THRESHOLD``
    Log2 fold change threshold for DE analysis (default: 0.1)

``--pval_threshold PVAL_THRESHOLD``
    P-value threshold for DE analysis (default: 0.05)

Enrichment Analysis Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--skip_enrichment``
    Skip enrichment analysis after limma (default: False)

``--enrich_databases ENRICH_DATABASES``
    Directory containing enrichment analysis databases (.gmt files) (default: species-specific, for example /omics_downstream_pipeline/ref_files/human/human_enrichr_databases/)

``--logfc_enrich LOGFC_ENRICH``
    Log2 fold change threshold for enrichment analysis. If not specified, ``--logfc_threshold`` value is used (default: value passed to logfc_threshold)

``--pval_enrich PVAL_ENRICH``
    P-value threshold for enrichment analysis. If not specified, ``--pval_threshold`` value is used (default: value passed to pval_threshold)

``--pval_enrich_column {adj.P.Val,P.Value}``
    P-value column to use for enrichment analysis either ``adj.P.Val`` or ``P.Value`` (default: "adj.P.Val")

``--min_size MIN_SIZE``
    Minimum allowed number of genes from gene set that are also in the dataset in gseapy.prerank function (default: 15)

``--max_size MAX_SIZE``
    Maximum allowed number of genes from gene set that are also in the data set in gseapy.prerank function (default: 1000)

``--permutation_num PERMUTATION_NUM``
    Number of permutations in gseapy.prerank function (default: 100)

Usage Examples
--------------

Complete Analysis with Batch Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    python DE_pipelines/TMT/tmt_command_line_script.py \
        --analysis_type both \
        --organism homo_sapiens \
        --expression_data_path data/protein_expression.csv \
        --protein_metadata_path data/protein_metadata.csv \
        --sample_metadata_path data/sample_metadata.csv \
        --preprocessing_output_dir results/preprocessing/ \
        --de_output_dir results/de_analysis/ \
        --condition_pairs control treatment \
        --batch_correct \
        --batch_id_column Batch \
        --mod_id_column Condition \
        --verbose

Preprocessing with Canonical Protein Mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    python DE_pipelines/TMT/tmt_command_line_script.py \
        --analysis_type preprocessing \
        --organism homo_sapiens \
        --expression_data_path data/protein_expression.csv \
        --protein_metadata_path data/protein_metadata.csv \
        --sample_metadata_path data/sample_metadata.csv \
        --preprocessing_output_dir results/preprocessing/ \
        --identify_canonical_proteins \
        --swiss_prot_ref_path ref_files/swiss_prot.txt \
        --apply_uniparc_mapping \
        --uniparc_db_path ref_files/uniparc.db

Differential Expression Only
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    python DE_pipelines/TMT/tmt_command_line_script.py \
        --analysis_type de \
        --input_file results/preprocessing/proteomics_data.h5ad \
        --de_output_dir results/de_analysis/ \
        --condition_pairs control treatment \
        --anndata_data_layer mad_log1p_sum_norm_peptide_norm \
        --ignore_non_canonical \
        --logfc_threshold 0.2 \
        --pval_threshold 0.01

Multiple Organisms Support
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # For mouse data
    python DE_pipelines/TMT/tmt_command_line_script.py \
        --analysis_type both \
        --organism mus_musculus \
        --expression_data_path data/mouse_protein_expression.csv \
        --protein_metadata_path data/mouse_protein_metadata.csv \
        --sample_metadata_path data/mouse_sample_metadata.csv \
        --preprocessing_output_dir results/mouse_preprocessing/ \
        --de_output_dir results/mouse_de_analysis/ \
        --condition_pairs wildtype knockout

Advanced Normalization Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Skip specific normalization steps
    python DE_pipelines/TMT/tmt_command_line_script.py \
        --analysis_type preprocessing \
        --organism homo_sapiens \
        --expression_data_path data/protein_expression.csv \
        --protein_metadata_path data/protein_metadata.csv \
        --sample_metadata_path data/sample_metadata.csv \
        --preprocessing_output_dir results/preprocessing/ \
        --skip_sum_normalize \
        --skip_peptide_normalize \
        --target_sum 500000 \
        --umap_plot \
        --n_neighbors 10

Custom Configuration with JSON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    python DE_pipelines/TMT/tmt_command_line_script.py \
        --analysis_type both \
        --preprocessing_config config/tmt_preprocessing_config.json \
        --de_config config/tmt_de_config.json \
        --verbose

Excel Output Control
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Handle large result sets that exceed Excel limits
    python DE_pipelines/TMT/tmt_command_line_script.py \
        --analysis_type both \
        --organism homo_sapiens \
        --expression_data_path data/protein_expression.csv \
        --protein_metadata_path data/protein_metadata.csv \
        --sample_metadata_path data/sample_metadata.csv \
        --preprocessing_output_dir results/preprocessing/ \
        --de_output_dir results/de_analysis/ \
        --condition_pairs control treatment \
        --ignore_excel_limit \
        --skip_plot_oxphos \
        --verbose