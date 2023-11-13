# N.B: You should create a separate directory for each run or experiment and
# copy/run this script and 'nygc.config' (which lets you run via slurm) in that
# directory.

# This is the path to the shared container images (NOT WRITABLE)
export NXF_SINGULARITY_LIBRARYDIR=/gpfs/commons/groups/imielinski_lab/data/pipeline/container_images_cache
# This is the path to your local container images (WRITABLE)
export NXF_SINGULARIY_CACHEDIR="${HOME}/local_singularity_cache"
mkdir -p $NXF_SINGULARIY_CACHEDIR

module unload java && module load java;
module load singularity/3.8.6;

mkdir hg19  # This is necessary for fragcounter to work

# Currently pointing to the 'mskilab' branch of the pipeline repo
# You should change this to your pipeline path if you want to run your own
# modified version of the pipeline
pipeline_dir=/gpfs/commons/groups/imielinski_lab/projects/nf-jabba

# You can specify the samplesheet.csv via the commandline
# Otherwise it will default to './samplesheet.csv'
samplesheet_csv=${1:-samplesheet.csv}

# Local MONSTER PON is specified to avoid downloading large file from AWS
# You can switch this to a different local PON if necessary
pon_path=/gpfs/commons/groups/imielinski_lab/data/dryclean/MONSTER_PON_RAW/MONSTER_PON_RAW_SORTED/fixed.detergent.rds

# You can change the flags/parameters below to suit your run.

# Use `nf-core launch` on the command line to launch a web based gui
# to modify any of the pipeline parameters from their defaults and then run

# Important parameters (see README.md for more details)
# --tools: tools listed are the minimal set necessary for JaBbA output
# --step: change this to 'alignment' if you are starting from fastqs
# -resume: always tries to resume if possible, otherwise will start from beginning
# -c nygc.config: include this flag/parameter if you want to run via slurm
nextflow run "$pipeline_dir" --input "$samplesheet_csv" \
	--outdir ./results/	\
	--genome GATK.GRCh37	\
	-profile singularity	\
	-with-report "report_$(date +'%Y%m%d_%H%M%S').html"	\
	-with-trace	\
	--tools svaba,hetpileups,fragcounter,dryclean,ascat,cbs,jabba	\
	--step sv_calling	\
	--pon_dryclean "$pon_path" \
    -c nygc.config \
	-resume

