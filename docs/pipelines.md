# Introduction to HCBC pipelines

Content - ![](https://img.shields.io/badge/status-WorkInProgress-yellow)

## Data Management

### AWS 

- Samplesheet input files for pipelines
    - `pipelineName_PI_hbcNNNNNN`
    - Have a copy in project folder in O2
    - Manually removing weekly during platform meeting
- Raw data is under `input` folder 
    - See instructions below to move data in/out
    - `pipelineName_PI_hbcNNNNNN`
    - lifecycle 14 days
- Pipeline outputs are under `results`:
    - `pipelineName_PI_hbcNNNNNN`
    - lifecycle 14 days for bigger than 1gb
    - Move output pipeline to project folder under `final` folder

### Move data in/out of AWS

Follow this to copy data in and out of our AWS space:

- Log in into transfer node in O2
- Type `sudo -su bcbio` to be bcbio user
- Use this command to copy data to AWS:
```
/usr/local/bin/aws s3 sync $FOLDER_WITH_FASTQ s3://hcbc-seqera/input/rnaseq_piname_hbcNNNN
```
- Use this command to copy data from AWS:
```
/usr/local/bin/aws s3 sync s3://hcbc-seqera/results/rnaseq_piname_hbcNNNN $FOLDER_PROJECT
```
**Make sure bcbio group has read/write access to the folders otherwise `aws` command won't work, but won't error either.**

## Parameters

### RNAseq

- We use salmon with bam files produced by STAR mapped to transcriptome for quantification

### CHIPseq

- It can analyze multiple antibodies in one pipeline run (pipeline splits samples by antibody)
- Default parameters
- de-duplication for all samples
- `bowtie` is set up with these extra parameters: `--sensitive-local -X 1000` (this is only true in seqera dev environment, not production)
- `macs_gsize` needs to be setup for each species accordingly tools

### CUT&RUN

- Run once per antibody (because pipeline does not split samples by antibody)
- Turn on `dedup_target_reads`
- Use both `macs2` and `seacr` for peakcalling (list macs2 first so it is used as primary)
- Normalization mode is set to `CPM` (can be changed if client has spike-in samples)
- Depending on the number of samples, user may want to skip `deeptools` processes involving all samples
- processes including SAMTOOLS_SORT, BEDTOOLS_SORT, SAMTOOLS_CUSTOMVIEW, FRAG_LEN_HIST, and DEEPTOOLS_PLOTHEATMAP_GENE_ALL are given more memory than nf-core default

### ATACseq

All peaks `nf-core-atac-seq_shift`:
- shift is on
- keep_dup is false

NFR peaks `nf-core-atac-seq_shift_NFR`:
- same than previous except parameters for Aligmentsieve:
  - `--minFragmentLength 0`
  - `--maxFragmentLength 120`

**Note**: Recommendation to check the fragment length distribution after the run to make sure you're capturing the NFRs.
**Note**: We don't need the MN(180, 247), DN (315, 473) and TN (558, 615), unless it's a specific case where we are looking at global shifts in accessibility

## Nextflow in Seqera platform

- Create an user here: https://cloud.seqera.io/login
- Ask Platform team to add you to HCBC workspace
- Transfer data to HCBC S3: Ask Alex/Lorena. Files will be at our S3 bucket `input/pipelineName_PI_hbcNNNNNN` folder

### RNAseq

- Prepare the CSV file according this [instructions](https://nf-co.re/rnaseq/3.14.0/docs/usage#multiple-runs-of-the-same-sample). File should look like this:

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,s3path/AEG588A1_S1_L002_R1_001.fastq.gz,s3path/AEG588A1_S1_L002_R2_001.fastq.gz,auto
CONTROL_REP1,s3path/AEG588A1_S1_L003_R1_001.fastq.gz,s3path/AEG588A1_S1_L003_R2_001.fastq.gz,auto
CONTROL_REP1,s3path/AEG588A1_S1_L004_R1_001.fastq.gz,s3path/AEG588A1_S1_L004_R2_001.fastq.gz,auto
```

Use `bcbio_nfcore_check(csv_file)` to check the file is correct.

You can add more columns to this file with more metadata, and use this file as the `coldata` file in the templates.

- Safe the file under `meta` folder
- Upload this file to our `Datasets` in Seqera using the name of the project but starting with `pipelineName_PI_hbcNNNNNN`
- Go to `Launchpad`, select `nf-core_rnaseq` pipeline, and select the previous created `Datasets` in the `input` parameter after clicking in `Browser`
  - Select an output directory with the same name used for the `Dataset` inside the `results/pipelineName_PI_hbcNNNNNN` folder in S3
- When pipeline is done, data will be copied to our on-premise HPC in the scratch system under `scratch/groups/hsph/hbc/bcbio/` folder


## Nextflow in O2

- Nextflow is available at `/n/app/bcbio/nextflow/nextflow`. 
- Singularity containers at available at `/n/app/singularity/containers/shared/bcbio/`.
- Cluster config: `/n/app/bcbio/nextflow/o2.config`

An example of sbatch script is:

```
#!/bin/bash

#SBATCH --job-name=Nextflow      # Job name
#SBATCH --partition=priority            # Partition name
#SBATCH --time=1-23:59                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=1            # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=12G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)

module load java/jdk-21.0.2
export NXF_APPTAINER_CACHEDIR=/n/app/singularity/containers/shared/bcbio/nf-core-rnaseq-3.14.0
export NXF_SINGULARITY_LIBRARYDIR=/n/app/singularity/containers/shared/bcbio/nf-core-rnaseq-3.14.0

/n/app/bcbio/nextflow/nextflow run nf-core/rnaseq -r 3.14.0 -profile singularity \
    -c /n/app/bcbio/nextflow/o2.config -c /n/app/bcbio/nextflow/rnaseq.resources.config \
    -params-file /n/app/bcbio/nextflow/rnaseq.json \
    --input samplesheet.csv --outdir this_folder -resume
```

### RNAseq

Containers at `/n/app/singularity/containers/shared/bcbio/nf-core-rnaseq-3.14.0`

### viralrecon

Read documentation [here](https://nf-co.re/viralrecon/2.6.0/). 


This is an example for test data:

```
module load java/jdk-21.0.2
export NXF_APPTAINER_CACHEDIR=/n/app/singularity/containers/shared/bcbio/nf-core-viralrecon_2.6.0
export NXF_SINGULARITY_LIBRARYDIR=/n/app/singularity/containers/shared/bcbio/nf-core-viralrecon_2.6.0

/n/app/bcbio/nextflow/nextflow run nf-core/viralrecon -r 2.6.0 -profile singularity,test --outdir this_folder -resume

```
To run your data, prepare input file following this [doc](https://nf-co.re/viralrecon/2.6.0/docs/usage/#samplesheet-format), and run it like this:

```
/n/app/bcbio/nextflow/nextflow run nf-core/viralrecon -r 2.6.0 -profile singularity --outdir this_folder --input samplesheet.csv -resume
```

## Nextflow in FAS

We use FAS to run pipelines in scratch, main storage location for data is still O2. 

- Run pipelines in scratch: `/n/netscratch/hsph_bioinfo/Lab`
- Keep downstream analysis in PIs folder: `/n/holylfs05/LABS/hsph_bioinfo/Lab/PIs`

Pipelines that have been run so far in FAS:
- rnaseq
- scrnaseq
- cutandrun

Before using nextflow, you need to load a recent version of java:

```
module load jdk/21.0.2-fasrc01
```

Use nextflow at `/n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/nextflow`

Use config file at `/n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/cannon.config`

Example command to run in an interactive job:
Note that this is only for test datasets requiring minimal CPUs, memory, and parallelization. Otherwise, submit as sbatch. 
```
/n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/nextflow run nf-core/rnaseq -profile test,singularity --outdir tmp -c /n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/fas.config
```

For non-test data, you will use sbatch to submit the head job of the pipeline, which will in turn submit the child jobs. 
Modify this template as needed before using it:

```
#!/bin/bash

#SBATCH --job-name=Nextflow      # Job name
#SBATCH --partition=shared            # Partition name
#SBATCH --time=0-48:59                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --mem=16G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)

module load jdk/21.0.2-fasrc01

# modify these paths as necessary to point to the containers for the pipeline you're using
export NXF_APPTAINER_CACHEDIR=/n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/nfcore-rnaseq
export NXF_SINGULARITY_LIBRARYDIR=/n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/nfcore-rnaseq

# Optional: if you'd like to monitor your run in Seqera Platform, set up a token there and use it here
export TOWER_WORKSPACE_ID=268530979103043
export TOWER_ACCESS_TOKEN=<your_access_token>

OUTPUT=path_to_results

/n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/nextflow run nf-core/rnaseq
  -r 3.14.0 \
  -profile singularity \
  # next line is for passing pipeline parameters. can also insead use --params-file and use JSON downloaded from seqera
  -c analysis.config \
  # next line is for optimizing resource requests, file does not yet exist for all pipelines
  -c rnaseq.resources.config \ 
  --outdir $OUTPUT
  -c cannon.config \
  # next line optional, if you want to monitor your run in Seqera Platform
  -with-tower
  -resume
```
