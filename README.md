# Notices

Git mirrors:
- [Codeberg](https://codeberg.org/paveloom-university/Computational-Workshop-S09-2021)
- [GitHub](https://github.com/paveloom-university/Computational-Workshop-S09-2021)
- [GitLab](https://gitlab.com/paveloom-g/university/s09-2021/computational-workshop)

Check the materials in each directory for assignments and other sources.

The reports are expected to be compiled with [`tectonic`](https://tectonic-typesetting.github.io/en-US/) as follows:

```bash
tectonic -X compile report.tex
```

This project provides [Julia](https://julialang.org) scripts. Make sure to use the project files (`Project.toml`) when running them:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/script.jl
```
