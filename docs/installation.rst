Installation
============

Requirements
------------
* Python 3.8+
* R 4.1+ (for DESeq2 and limma)
* Conda or Mamba package manager (recommended)

Installation Methods
-------------------

Method 1: Using Conda/Mamba (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Create conda environment**

   .. code-block:: bash

      mamba create -n name_of_environment -c conda-forge -c bioconda python=3.8 matplotlib seaborn "pandas<2.0.1" numpy anndata scipy statsmodels scikit-learn plumbum samtools star bedtools r-base=4.1 abra2 art ucsc-bedGraphToBigWig bowtie2 subread gatk4 piper=0.12 pysam trimmomatic rpy2 tzlocal ReportLab pytest-cov codecov libxml2 ipykernel jupyter gseapy fastcluster 'scanpy>1.9.3' r-xml r-mass r-matrix  bioconductor-tximport bioconductor-enhancedvolcano plotly dash dash-bio bioconductor-limma bioconductor-deseq2 jupyter ipykernel

2. **Activate the environment**

   .. code-block:: bash

      conda activate name_of_environment

3. **Clone the repository**

   .. code-block:: bash

      git clone https://github.com/MoothaLab/omics_downstream_pipeline.git
      cd omics_downstream_pipeline

4. **Install additional dependencies**

   .. code-block:: bash

      git clone git@github.com:wckdouglas/diffexpr.git
      cd diffexpr
      python setup.py install
      cd ..

Method 2: Using Conda Environment File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Clone the repository**

   .. code-block:: bash

      git clone https://github.com/MoothaLab/omics_downstream_pipeline.git
      cd omics_downstream_pipeline

2. **Create environment from file**

   .. code-block:: bash

      conda env create -f conda_environment.yml
      conda activate [environment_name]

3. **Install additional dependencies**

   .. code-block:: bash

      git clone git@github.com:wckdouglas/diffexpr.git
      cd diffexpr
      python setup.py install
      cd ..

Required Dependencies
--------------------

Core Packages
~~~~~~~~~~~~~
* **Python packages:** pandas, numpy, scipy, matplotlib, seaborn
* **Statistical analysis:** scanpy, anndata, statsmodels, scikit-learn
* **Bioinformatics:** pysam, gseapy, fastcluster

R Packages
~~~~~~~~~~
* **Differential expression:** DESeq2, limma
* **Visualization:** EnhancedVolcano
* **Data import:** tximport
* **Base packages:** XML, MASS, Matrix

Bioinformatics Tools
~~~~~~~~~~~~~~~~~~~
* **Alignment:** STAR, bowtie2
* **Processing:** samtools, bedtools, subread
* **Variant calling:** GATK4, abra2
* **Utilities:** UCSC tools, trimmomatic

Repository Structure
-------------------

.. code-block:: text

   omics_downstream_pipeline/
   ├── DE_pipelines/
   │   ├── RNAseq/              # RNA-seq differential expression analysis
   │   │   └── run_pipeline_command_line.py
   │   ├── TMT/                 # TMT proteomics analysis  
   │   │   └── tmt_command_line_script.py
   │   └── utilities/           # Shared utility functions
   ├── docs/                    # Documentation and images
   │   └── pipeline_overview.png
   ├── ref_files/               # Reference files and databases
   ├── conda_environment.yml   # Conda environment specification
   └── README.md               # Project documentation

Troubleshooting
--------------

Environment Issues
~~~~~~~~~~~~~~~~~
If you encounter issues with package conflicts, try creating a fresh environment:

.. code-block:: bash

   conda deactivate
   conda env remove -n name_of_environment
   # Then repeat Method 1 installation steps

R Package Issues
~~~~~~~~~~~~~~~
If R packages fail to install via conda, you can install them manually in R:

.. code-block:: r

   if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   BiocManager::install(c("DESeq2", "limma", "EnhancedVolcano", "tximport"))

Verification
-----------
To verify your installation, test importing key packages:

.. code-block:: python

   import pandas as pd
   import numpy as np
   import scanpy as sc
   import anndata as ad
   import gseapy as gp
   
   # Test R integration
   import rpy2.robjects as robjects
   from rpy2.robjects import pandas2ri
   
   print("Installation successful!")

Repository Structure
-------------------

.. code-block:: text

   omics_downstream_pipeline/
   ├── DE_pipelines/
   │   ├── RNAseq/              # RNA-seq differential expression analysis
   │   │   └── run_pipeline_command_line.py
   │   ├── TMT/                 # TMT proteomics analysis  
   │   │   └── tmt_command_line_script.py
   │   └── utilities/           # Shared utility functions
   ├── docs/                    # Documentation and images
   │   └── pipeline_overview.png
   ├── ref_files/               # Reference files and databases
   └── README.md               # This documentation
