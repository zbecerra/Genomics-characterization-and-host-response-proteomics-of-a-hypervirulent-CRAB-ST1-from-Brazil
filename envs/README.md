# Conda environments

The tools were run in separate conda environments (a tool → environment table is in [`docs/requirements.md`](../docs/requirements.md)).
Generate the environment files once, on the cluster where they were used, and commit them here:

```bash
bash envs/export_envs.sh          # writes envs/<name>.yml for every environment found
```

Recreate one environment elsewhere with `conda env create -f envs/<name>.yml`.
Notes: `mlst` needs a patch after installation (`scripts/00_setup/patch_mlst_blast_check.sh`); the genome-wide AbaR island screen and the ISfinder search
used BLAST+ 2.5.0+, the flanking-region comparison BLAST+ 2.16.0 (see `docs/requirements.md`).

Environment files were exported with `conda env export --no-builds`; machine-specific `prefix:` lines were removed.
Two conda installations were used on the cluster (`~/.conda/envs` and `~/miniconda3/envs`), and Bakta lives in a third, system-wide one; `bakta.yml` is therefore provided as a package list (`bakta_conda_list.txt`) once exported (see the repository checklist).
