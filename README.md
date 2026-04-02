
# Claude-Mendelian-Randomization

Developed on Claude


Download select reference files...

```bash

./download_icvf_sumstats.sh

```

Install needed R libraries ...

```bash
module load r

Rscript -e 'install.packages("data.table")'
Rscript -e 'remotes::install_github("MRCIEU/TwoSampleMR")'
Rscript -e 'remotes::install_github("MRCIEU/genetics.binaRies")'

```

Format the references files for processing ...

```bash

sbatch --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 --export=None \
  --job-name=format --ntasks=1 --cpus-per-task=64 --mem=490G  --output="format_icvf_for_mr.log" \
  --wrap="module load r; format_icvf_for_mr.R"

```

Get your own token file from https://api.opengwas.io/ to run the next step.

Run MR ...

```bash

sbatch --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 --export=None \
  --job-name=mr --ntasks=1 --cpus-per-task=64 --mem=490G  --output="run_bidirectional_mr.log" \
  --wrap="module load r plink; run_bidirectional_mr.R"

```


Download the pre-built MRCIEU reference panel — this is literally the same data the IEU API uses internally, already in plink format:

```bash

wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
tar -xzf 1kg.v3.tgz

```

Ensure that this path is correct in the following script.

Post process ...

```bash

sbatch --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 --export=None \
  --job-name=post --ntasks=1 --cpus-per-task=64 --mem=490G  --output="mr_postprocessing.log" \
  --wrap="module load r plink; mr_postprocessing.R"

```



