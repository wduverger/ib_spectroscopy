# Single-cell infrared absorption spectroscopy of bacterial inclusion bodies

[![CC BY 4.0][cc-by-shield]][cc-by-link]
[![DOI (repository)][doi-repo-shield]][doi-repo-link]
[![DOI (paper)][doi-paper-shield]][doi-paper-link]

![Abstract Graphic](./figures/abstract%20graphic.png)

This repository contains the data and analysis presented in **Duverger et al. 2024. (submitted)**

The analysis notebooks can be found in the [./code]() directory. They provide all figures and results presented in the manuscript. The underlying data can be found in the [./primary_data]() directory.

The code can be ran by executing the following commands, assuming git and conda are installed.

```
    git clone https://github.com/wduverger/ib_spectroscopy
    cd ib_spectroscopy/code
    conda env create -f ./environment.yml -p ./.conda_env
    conda activate ./.conda_env
    jupyter lab
```

## Citing this work

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by-link].
If you use this code in your work, please cite the following publication:

> Duverger et al. Single-cell infrared absorption spectroscopy of bacterial inclusion bodies. 2024. (submitted)

When explicitly referencing the data or code, use the DOI assigned to this repository: [DOI forthcoming][doi-repo-link]

[cc-by-link]:       http://creativecommons.org/licenses/by/4.0/
[cc-by-image]:      https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]:     https://img.shields.io/badge/License-CC%20BY%204.0-tomato.svg
[doi-paper-shield]: https://img.shields.io/badge/DOI%20(paper)-forthcoming-tomato
[doi-repo-shield]:  https://img.shields.io/badge/DOI%20(repository)-forthcoming-tomato
[doi-paper-link]:   https://doi.org/forthcoming
[doi-repo-link]:    https://doi.org/forthcoming