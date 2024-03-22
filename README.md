# Single-cell infrared absorption spectroscopy of bacterial inclusion bodies

[![CC BY 4.0][cc-by-shield]][cc-by-link]
[![DOI (repository)][doi-repo-shield]][doi-repo-link]
[![DOI (preprint)][doi-paper-shield]][doi-paper-link]

![Abstract Graphic](./figures/abstract%20graphic.png)

This repository contains the data and analysis presented in **Duverger et al. 2024. (submitted)**

The analysis notebooks can be found in the [code](./code) directory. They provide all figures and results presented in the manuscript. The underlying data can be found in [primary_data](./primary_data).

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
[Creative Commons Attribution 4.0 International License (CC BY 4.0)][cc-by-link].
If you use this code in your work, please cite the following preprint:

> Wouter Duverger, Grigoria Tsaka, Ladan Khodaparast et al. Single-cell infrared absorption spectroscopy of bacterial inclusion bodies, 19 March 2024, PREPRINT (Version 1) available at Research Square [10.21203/rs.3.rs-4094624/v1][doi-paper-link]

When explicitly referencing the data or code, use the DOI assigned to this repository: [10.6084/m9.figshare.25398622][doi-repo-link]

[cc-by-link]:       http://creativecommons.org/licenses/by/4.0/
[cc-by-image]:      https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]:     https://img.shields.io/badge/License-CC%20BY%204.0-tomato.svg

[doi-paper-shield]: https://img.shields.io/badge/DOI%20(preprint)-10.21203/rs.3.rs--4094624/v1-tomato
[doi-paper-link]:   https://doi.org/10.21203/rs.3.rs-4094624/v1

[doi-repo-shield]:  https://img.shields.io/badge/DOI%20(repository)-10.6084/m9.figshare.25398622-tomato
[doi-repo-link]:    https://doi.org/10.6084/m9.figshare.25398622
