Usage Overview
==============

Quick Start
-----------

Both RNA-seq and TMT pipelines support:

* Command-line execution with extensive parameters
* JSON configuration files for reproducible analyses
* Modular execution (preprocessing only, DE only, or both)
* Multi-organism support (homo_sapiens, mus_musculus, other)
* Advanced normalization and batch correction options

Basic Workflow
--------------

1. **Set up environment**: Install dependencies using conda/mamba
2. **Prepare your data**: Organize count files and metadata according to pipeline requirements
3. **Run preprocessing**: Quality control, normalization, and visualization
4. **Differential expression**: Statistical analysis between conditions using DESeq2 (RNA-seq) or Limma (TMT)
5. **Enrichment analysis**: Functional annotation of results using pathway databases

Pipeline Selection Guide
-----------------------

RNA-seq Pipeline
~~~~~~~~~~~~~~~

Use for:
* Gene expression analysis from RNA sequencing
* featureCounts or similar count data
* Large-scale transcriptomic studies
* Time-series or multi-condition experiments

Key features:
* DESeq2-based differential expression
* Support for protein_coding and lncRNA transcripts
* OXPHOS complex analysis
* PCA and UMAP visualization

TMT Proteomics Pipeline
~~~~~~~~~~~~~~~~~~~~~~

Use for:
* Tandem Mass Tag (TMT) proteomics data
* Protein abundance measurements
* Batch effect correction capabilities
* Canonical protein isoform identification

Key features:
* Limma-based differential expression
* Multiple normalization strategies (sum, peptide, MAD)
* Combat batch correction
* UniParc and Swiss-Prot integration

Data Requirements
----------------

RNA-seq Data Structure
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   project/
   ├── data/
   │   ├── sample1/
   │   │   └── featureCounts.txt
   │   ├── sample2/
   │   │   └── featureCounts.txt
   │   └── ...
   ├── metadata/
   │   └── sample_metadata.csv
   └── results/

**Required metadata columns:**
* ``RNA-seq Ref``: Sample names (string)
* ``Replicate``: Replicate information (string)
* ``Condition``: Experimental conditions
* Additional columns for covariates, batch info, etc.

TMT Data Structure
~~~~~~~~~~~~~~~~~

.. code-block:: text

   project/
   ├── data/
   │   ├── protein_expression.csv
   │   ├── protein_metadata.csv
   │   └── sample_metadata.csv
   └── results/

**Required metadata columns:**
* ``TMT ID``: Sample names (string)
* ``Replicate``: Replicate information (string)
* ``Condition``: Experimental conditions
* Additional columns for batch correction, covariates, etc.

Configuration Options
--------------------

Command Line Parameters
~~~~~~~~~~~~~~~~~~~~~~

Both pipelines accept extensive command-line parameters for customization:

.. code-block:: bash

   # View all available options
   python DE_pipelines/RNAseq/run_pipeline_command_line.py -h
   python DE_pipelines/TMT/tmt_command_line_script.py -h

JSON Configuration Files
~~~~~~~~~~~~~~~~~~~~~~~~

For reproducible analyses, use JSON configuration files:

.. code-block:: bash

   # RNA-seq with configuration files
   python DE_pipelines/RNAseq/run_pipeline_command_line.py \
       --analysis_type both \
       --preprocessing_config config/rna_preprocessing.json \
       --de_config config/rna_de.json

   # TMT with configuration files
   python DE_pipelines/TMT/tmt_command_line_script.py \
       --analysis_type both \
       --preprocessing_config config/tmt_preprocessing.json \
       --de_config config/tmt_de.json

Organism Support
~~~~~~~~~~~~~~~

Both pipelines support multiple organisms:

* ``homo_sapiens``: Human (default reference files)
* ``mus_musculus``: Mouse
* ``other``: Custom organism (requires manual reference file specification)

.. code-block:: bash

   # Example for mouse data
   python DE_pipelines/RNAseq/run_pipeline_command_line.py \
       --analysis_type both \
       --organism mus_musculus \
       --feature_counts_pattern "mouse_data/*/featureCounts.txt" \
       --sample_metadata_path metadata/mouse_samples.csv

Analysis Types
-------------

Preprocessing Only
~~~~~~~~~~~~~~~~~

Run quality control, normalization, and visualization:

.. code-block:: bash

   # RNA-seq preprocessing
   python DE_pipelines/RNAseq/run_pipeline_command_line.py \
       --analysis_type preprocessing \
       --organism homo_sapiens \
       --feature_counts_pattern "data/*/featureCounts.txt" \
       --sample_metadata_path metadata/samples.csv \
       --preprocessing_output_dir results/preprocessing/

   # TMT preprocessing with batch correction
   python DE_pipelines/TMT/tmt_command_line_script.py \
       --analysis_type preprocessing \
       --organism homo_sapiens \
       --expression_data_path data/protein_expression.csv \
       --protein_metadata_path data/protein_metadata.csv \
       --sample_metadata_path data/sample_metadata.csv \
       --batch_correct \
       --batch_id_column Batch

Differential Expression Only
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run statistical analysis on preprocessed data:

.. code-block:: bash

   # RNA-seq DE analysis
   python DE_pipelines/RNAseq/run_pipeline_command_line.py \
       --analysis_type de \
       --input_file results/preprocessing/rna-seq_data.h5ad \
       --output_dir results/de_analysis/ \
       --condition_pairs control treatment

   # TMT DE analysis
   python DE_pipelines/TMT/tmt_command_line_script.py \
       --analysis_type de \
       --input_file results/preprocessing/proteomics_data.h5ad \
       --de_output_dir results/de_analysis/ \
       --condition_pairs control treatment

Complete Pipeline
~~~~~~~~~~~~~~~~

Run both preprocessing and differential expression:

.. code-block:: bash

   # Complete RNA-seq analysis
   python DE_pipelines/RNAseq/run_pipeline_command_line.py \
       --analysis_type both \
       --organism homo_sapiens \
       --feature_counts_pattern "data/*/featureCounts.txt" \
       --sample_metadata_path metadata/samples.csv \
       --preprocessing_output_dir results/preprocessing/ \
       --output_dir results/de_analysis/ \
       --condition_pairs control treatment \
       --verbose

Output Files
-----------

Both pipelines generate:

* **AnnData files** (``.h5ad``): Processed data for downstream analysis
* **Statistical results**: Differential expression tables
* **Visualizations**: PCA plots, heatmaps, volcano plots
* **Quality control**: Summary statistics and diagnostic plots
* **Enrichment results**: Pathway analysis (if not skipped)

Large Dataset Handling
---------------------

For datasets with many results that exceed Excel display limits:

.. code-block:: bash

   # Use ignore_excel_limit flag
   python DE_pipelines/RNAseq/run_pipeline_command_line.py \
       --analysis_type both \
       --ignore_excel_limit \
       --condition_pairs control treatment

   python DE_pipelines/TMT/tmt_command_line_script.py \
       --analysis_type both \
       --ignore_excel_limit \
       --condition_pairs control treatment

Command Line Tips
-----------------

Essential Flags
~~~~~~~~~~~~~~

* Always run scripts with ``-h`` flag first to see all available options
* Use ``--verbose`` flag for detailed output during analysis
* Use ``--organism`` to specify your species for appropriate reference files
* Use ``--ignore_excel_limit`` for large result sets

File Organization
~~~~~~~~~~~~~~~~

* Reference files are automatically loaded from ``ref_files/`` directory
* Output files are saved in AnnData (.h5ad) format for interoperability
* Organize input data according to pipeline-specific requirements
* Use consistent naming conventions for samples and conditions

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~

* Use ``--skip_plot_oxphos`` to speed up analysis if OXPHOS analysis not needed
* Consider using configuration files for complex parameter sets
* Run preprocessing and DE separately for large datasets to checkpoint progress

Troubleshooting
--------------

Common Issues
~~~~~~~~~~~~

1. **Missing reference files**: Ensure organism-specific reference files exist
2. **Metadata column names**: Check required column names match exactly
3. **File format issues**: Verify CSV files are properly formatted
4. **Memory issues**: Consider running steps separately for very large datasets

Getting Help
~~~~~~~~~~~

* Check parameter documentation with ``-h`` flag
* Verify input file formats match requirements
* Use ``--verbose`` for detailed error messages
* Ensure all required dependencies are installed

