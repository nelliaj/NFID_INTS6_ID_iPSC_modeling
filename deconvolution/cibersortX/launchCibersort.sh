#!/bin/bash -l
#SBATCH --account=project_deconvolution
#SBATCH --job-name=cibersortX
#SBATCH --output=log/csX_%A_%a.out
#SBATCH --error=log/csX_%A_%a.err
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=60G


INPUTLINE=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n1)
echo $INPUTLINE
scRNAseqRef=$(echo $INPUTLINE | awk {'print $1'})
bulkQuery=$(echo $INPUTLINE | awk {'print $2'})
outputName=$(echo $INPUTLINE | awk {'print $3'})

apptainer exec -B deconvolution/cibersortx/input/:/src/data -B deconvolution/cibersortx/output/${outputName}/:/src/outdir deconvolution/apptainer/cibersortx-fractions.sif /src/CIBERSORTxFractions --username <registered_mail> --token <personalaccesstoken> --single_cell TRUE --refsample ${scRNAseqRef} --mixture ${bulkQuery} --fraction 0 --rmbatchSmode TRUE 

