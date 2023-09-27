# powr-refactor

![License: GPL-3](https://img.shields.io/github/license/ssciwr/powr-refactor)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/ssciwr/powr-refactor/ci.yml?branch=main)
![Language](https://img.shields.io/github/languages/top/ssciwr/powr-refactor)

This is a development repository for refactoring and improving the numerical stability of the coli part of the [PoWR code](https://github.com/powr-code/PoWR) (Potsdam Wolf-Rayet Stellar Atmospheres). For a description of PoWR and the available models, see [here](https://www.astro.physik.uni-potsdam.de/~wrh/PoWR/powrgrid1.php).

*This library is currently under development!*

## Workflow

The source code is a collection of >400 Fortran77 files with >500 subroutines. For a complete PoWR cycle, different programs are called consecutively. The execution is handled by bash scripts.

In this repository, the development focuses on the coli part of the execution cycle, in particular the [colimo](src/colimo.f) subroutine.

Testing needs to include one and several coli cycles, as well as integration into the PoWR execution cycle.

It also includes compilation of the source code on MacOS and Ubuntu OS with different intel and gnu compilers. The results from different compilers and different optimizations need to be consistent.

## Tasks
- initial set-up of the CI to ensure valid output
- improve the make process: better structure of source directory and make process
- improve the make process: include gnu compiler and debug compiler options
- profile the test runs for time and memory consumption
- implement Fortran best practices in `colimo.f` and modernize to latest standards
- identify computational and numerical bottlenecks and improve the algorithm


## Detailed description of the tests
To run all the commands, you need to source `.powrconfig` and `powr/proc.dir/bash_setup` to set the environment. 

The testmodel files can be downloaded using
```
wget -O testmodels.tgz https://heibox.uni-heidelberg.de/f/a62c7ae5559d43a0a8b2/?dl=1
```

### `coli` test
Create a new chain, ie 1, by `makechain 1`. This copies some scripts and executables into different folders:

| Folder      | Purpose |
| ----------- | ----------- |
| wrdata1     | data directory of the current results |
| scratch subdirectories | directories to handle execution cycle process |
| wrjobs | collection of scripts to run |
| output | directory containing the global output |
| tmp_data | intermediate results? commonly scratch? |
| tmp_2day | ? |

The most important files are in the `wrdata1` directory: `CARDS, DATOM, FEDAT, FEDAT_FORMAL, FGRID, FORMAL_CARDS, MODEL, MODEL_STUDY, NEWDATOM_INPUT, NEWFORMAL_CARDS_INPUT, next_job, next_jobz`  
Most of these are input files, some control the job execution.

The test is then run by `sub colitest1`. This creates the run log in `output` (`colitest1.log` and `colitest1.cpr`) - these are checked if the run are successful (`COLITEST  finished` in `log`). The results in `cpr` are compared to the reference file.

The output that is generated in `wrdata1` is also compared: `MODEL_STUDY_DONE, MODEL_STUDY_DONE_STEAL`


### Integration test

### Unit tests
TBD
