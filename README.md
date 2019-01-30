## A Large Scale, Cloud Based, Low Cost and Reproducible Mutation Calling Pipeline Using Docker Containers.

## Abstract.

The identification of variants that occur in the human genome remain a critical step in the analysis of Next Generation Sequencing (NGS) data. Accurate and timely variant calling is critical to precision medicine which seeks to relate particular genomic variations in a patient’s genome to targetable genes, drugs and treatments tailored to each patient. 

However, the analysis of these low-frequency variants require very sensitive algorithms within computational intensive pipelines. With high costs and luck of technology, most institutions resort to running bioinformatics pipelines as bash scripts on local clusters which are not only slow and cumbersome but also difficult to implement. This has in part led to reduce reproducibility and exhaustion of local data storage at such institutions. 

Therefore, I demonstrate how to implement bioinformatics pipelines in the cloud using open source tools like Github, Docker and Broad Institute Genome Analysis Toolkit (GATK). This fast and low cost implementation leverages parallel execution and auto scaling within the cloud to handle the computationally intensive pipeline. 

I show that such pipelines produce accurate results by analysing data from the 1000 Genomes Project. This would allow researchers to use more of their time doing research and less time configuring workflows.
