# isr_hsc_kinetics

Starter repository for a bioinformatics project with Claude Code running in a devcontainer. The biological question, datasets, and software license are pending.

## Start

1. Start Docker Desktop. Install VS Code and its Dev Containers extension.
2. Open this repository folder and run **Dev Containers: Reopen in Container**.
3. Wait for the environment setup and smoke check to finish.
4. In the container terminal, run `claude` and complete sign-in. Claude state persists in a project-specific Docker volume.
5. Select `.venv/bin/python` for Python scripts and notebooks. Run `bash scripts/check-environment.sh` to recheck the environment.

The container includes Python 3.12, JupyterLab, scientific Python packages, Biopython, Ubuntu R, samtools, bcftools, bedtools, tabix, FastQC, Git LFS, and GitHub CLI. Add Bioconductor and workflow-specific tools after selecting the analysis. No dataset is downloaded automatically.

## Layout

- `metadata/`: tracked dataset and sample provenance.
- `data/`: ignored raw, interim, and processed datasets.
- `scripts/`: reusable download, QC, and analysis commands.
- `workflows/`: future workflow definitions.
- `notebooks/`: exploratory analyses.
- `results/`: ignored generated outputs.

## Reproducibility

This is a bootstrap environment, not a frozen scientific environment. Features, apt packages, and Python version ranges can resolve differently on rebuild. Before analysis, resolve Python dependencies into a committed lock file, initialize `renv` for R work, and record container image digests and tool versions. Set dataset-specific reference genomes, seeds, and parameters explicitly. The metadata manifest starts empty.

## Container boundaries

The project folder is writable from the container. The configuration does not mount the Docker socket, host home, SSH directory, or cloud credentials. VS Code may forward Git credentials or an SSH agent according to its own settings; disable forwarding if strict credential isolation is needed. Network access is enabled for public datasets and Claude. This is a development container, not a hardened sandbox; use normal Claude permission prompts.

## Keeping copies in sync

GitHub is the canonical remote for committed code. The devcontainer bind-mounts this local checkout, so host and container edit the same files without a separate sync step. The named Docker volume persists Claude state, not a second code checkout.

Before starting work, run `git pull --ff-only`. After reviewing changes, commit and run `git push`. Git does not automatically sync uncommitted edits. Ignored datasets, results, and Claude state are not backed up to GitHub; keep separate storage or reproducible download commands for those files.

## GitHub initialization

On the host, authenticate with `gh auth login`. After choosing a repository name and checking the files:

```bash
git add .
git commit -m "Initialize bioinformatics devcontainer project"
gh repo create OWNER/isr_hsc_kinetics --public --source=. --remote=origin --push
```

Choose an appropriate code license before inviting reuse; public visibility alone does not grant a software license. Dataset terms remain separate.

## References

- [Claude Code devcontainer documentation](https://code.claude.com/docs/en/devcontainer)
- [Development Container specification](https://containers.dev/)
