# Evaluation of funannotate

Analysed by: Lukas Jelonek

## Research questions/topics

1. Gain knowledge about the annotation process with `funannotate`
2. Evaluate if the process can be automized


## Research plan

1. Gather data for evaluation (contigs/scaffolds + rna-seq data)
2. Install `funannotate`
3. Annotate downloaded data with `funannotate`

## Data and Software

### P. cubensis

* Reasons: Personal interest, moderate size genome, rna-seq data available

* Assembly source: https://mycocosm.jgi.doe.gov/Psicub1_1/Psicub1_1.home.html
* RNA-Seq source: https://www.ebi.ac.uk/ena/browser/view/PRJNA450675

### funannotate

* Source: conda
* Documentation: https://funannotate.readthedocs.io/en/latest/

## Analysis

The whole analysis is described by the nextflow workflow in `main.nf`. This
workflow has a few dependencies that are not (yet) included and thus must be
installed/configured beforehands.

### Setup of funannotate

```bash
# install everything via conda
conda create --prefix <your-installation-target-directory> funannotate
conda activate --prefix <your-installation-target-directory>

# setup databases
export FUNANNOTATE_DB= /vol/funant/evaluate_funannotate/funannotate_dbs/

# install database by database as the interpro installation did not work
# and must be excluded
for i in merops uniprot dbCAN pfam repeats go mibig busco_outgroups gene2product; do funannotate setup -w -i $i; done
funannotate setup -i busco -b dikarya
```

### Setup of eggnogmapper

```bash
# funannotate had no option to configure the eggnog-mapper database-directory. I added 
# the possibility to configure the eggnog-mapper database-directory via an env-variable.
# It now has been merged to eggnog-mapper and should be available in future releases (
# see https://github.com/eggnogdb/eggnog-mapper/pull/220 )
export EGGNOG_DATA_DIR=/vol/funant/evaluate_funannotate/eggnog_db/
download_eggnog_data.py --data_dir $EGGNOG_DATA_DIR 
```

### Run the workflow

```bash
nextflow run . -resume -with-report report.html -with-dag graph.html
```

### Results

* `results/` - the funannotate annoation output
* `graph.html` - the workflow graph
* `report.html` - the runtime report of the nextflow workflow
