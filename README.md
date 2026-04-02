
# Claude-Operon-Mendelian-Randomization

Developed on Claude Operon



```bash

./download_icvf_sumstats.sh

```



```bash
module load r

Rscript -e 'install.packages("data.table")'
Rscript -e 'remotes::install_github("MRCIEU/TwoSampleMR")'

```


```bash

sbatch --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 --export=None --job-name=format --ntasks=1 --cpus-per-task=64 --mem=490G  --output="format.log" --wrap="module load r; format_icvf_for_mr.R"

```

You will need to get your own token file from https://api.opengwas.io/ to run the next step.

```bash

sbatch --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 --export=None --job-name=mr --ntasks=1 --cpus-per-task=64 --mem=490G  --output="mr.log" --wrap="module load r; run_bidirectional_mr.R"

```



