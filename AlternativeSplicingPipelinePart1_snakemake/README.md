# Sicelore2.0

## Installation

To install sicelore2.0 environnement
```
conda env create -n use_mamba mamba
conda activate use_mamba
mamba env create -f env.yml
```

If you don't have sicelore installed, then : 
```git clone https://github.com/ucagenomix/sicelore.git ```

## Getting started

#### Defining a snakemake profile
You can modify the snakemake profile, by changing the default resources in the cluster_config.yml file (see Snakemake documentation). Using this profile will create automatically a log folder with the logs of each individual rule.
Snakemake needs the profile to be stored at a specific place, so let's move this profile at the right location.

```
mkdir ~/.config/snakemake
mv profiles/profile_snakemake ~/.config/snakemake -R
```
#### Launching the pipeline
The only file that needs to be modified is the config.yml file.
If you want to use the default configuration (that works fine for ~200,000,000 reads in the original ONT fastq file), just modify the following:
- software paths
- sample_paths
- outputs
- num_split : you don't want to have fastq files with too few reads
- junction files : used for minimap in splice mode as a guidance for expected junctions
- raw_counts : with your own dataset. Obtaining this csv file can be done with the getting_raw_counts.R file in the utils folder

```
getting_raw_counts.R utilization :

Rscript getting_raw_counts.R path/to/output/cellranger/filtered_feature_bc_matrix/ path/to/mydir/raw_counts.csv
```

#### Expected inputs
```
|-- sample_folder
    |-- input_files
        |-- 1.ONT_fastq
            |-- ont.fastq

        |-- 2.short_read_files
            |-- possorted_genome.bam
            |-- possorted_genome.bam.bai
```

#### Outputs
Here's the structure of the output's directory, as long with the important output files.

```
|-- outputs_folder

    |-- sicelore_outputs_sample_name
    
        |-- consensus_outputs
        
            |-- consensus.sorted.tags.GE.bam 
                (original consensus bamfile)
                
            |-- consensus.final_output.sorted.tags.GE.bam 
                (added saved reads to the consensus computation)
                
        |-- GEUS10xAttributes.umifound.bam 
            (original output of sicelore) 
            
        |-- GEUS10xAttributes.all.bam 
            (bamfile with all the reads that map to the genome, no error correction)
            
        |-- GEUS10xAttributes.final_output.bam 
            (output of sicelore + saved reads)
            
        |-- isoform_outputs
            see sicelore pipeline
```
