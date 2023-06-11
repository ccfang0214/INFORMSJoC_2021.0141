[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Data and Codes for "A Study on Optimal Release Schedule for Multi-Version Software"

This repository includes the R programming codes for evaluating the Cost and Goodness-of-fit of Multi-Version SRGMs for the following paper:

Huang YS*, Fang CC, Chou CH, Tseng ZL, A Study on Optimal Release Schedule for Multi-Version Software. INFORMS Journal on Computing, 2023.

## Cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

To cite this software, please cite the paper using its DOI and the software itself, using the following DOI.

```
@article{DBOTSROA21,
  author =        {Y.S. Huang*, C.C. Fang, C.H. Chou, and Z.L. Tseng},
  publisher =     {INFORMS Journal on Computing},
  title =         {A Study on Optimal Release Schedule for Multi-Version Software},
  year =          {2023},
  doi =           {10.5281/zenodo.XXXXXXX},
  url =           {https://github.com/INFORMSJoC/2021.0141},
}  
```

## Description
The purpose of this repository is to share the data and codes for measuring the performance of goodness-of-fit analysis of multi-version software reliability growth models and evaluating the software testing cost under different cases. Our motivation is to present our application so that it can be easily replicated or further research can be conducted based on the analysis presented.

## Repository Structure
There are two folders for readers’ reference. These codes were written in R programming language. The folder “GoodnessofFit” includes seven code files for obtaining the least squares estimation of the related model parameters under different testing datasets. The figures of the fitting results are also generated from these program files. The folder “CostEvaluation” includes three code files for evaluating the software testing cost under different cases. The settings of related data and parameters are also included in these program files. 