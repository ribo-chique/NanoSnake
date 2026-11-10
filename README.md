# NanoSnake
Snakemake workflow for Nanopore direct RNA sequencing modification analysis. Made for small number of reference transcripts with well annotated modification sites.

## Description
A complete revamp of the sequencing workflow outlined [here](https://github.com/Thomson-RNA-Lab/Integrated-rRNA-snoRNA-seq)! Principle remain the same, by using IVT control sequences, we set a target sensitivity threshold in order to obtain a more stringent modification threshold for modkit pileup.
This is built for HPC with SLURM scheduler, but can be adapted without job scheduling.

## Key changes
- Automated run: Scripts for different steps (eg basecall, align) are all joined together. No more starting and stopping.
- Parallel job submission: Multiple samples are analysed together.
- Auto-processing of IVT samples: The script recognises IVT samples by sample name; multiple IVT replicates are detected and merged to obtain an average sample for thresholding to occur.
- More customisable parameters: Some steps, particularly for QC and modkit allow user-input parameters to supplement modification basecalling.

## Usage
Since this is built for HPC with SLURM scheduler, the set up for scheduler is optional
1. Create apptainer image with dorado and modkit
   ```apptainer build nanopore-essentials_copy.def```
2. Install conda env in `nanopore_smk.yml`
   ```conda env create -f nanopore_smk.yml```
3. `profile.yaml` stores information about default slurm scheduler limits, which is dependent of specific HPCs. Copy `profile.yaml` into `~/.config/snakemake/slurm`
4. Adjust config settings in `workflow/config/config.yaml`
5. Submit job!
   ```sbatch submit_snakemake```
 
## To be updated
There are more improvements I want to make but I don't have the time
- [ ] Select aligner in config: minimap2 & bwa
- [ ] Select modkit versions: v0.6 is a bit of a pain for high sequencing depth datasets. I would rather run v0.5 than to mess with sub-sampling `.bam` to have lower depth, pileup on sub-samples, then merge pileup outputs
