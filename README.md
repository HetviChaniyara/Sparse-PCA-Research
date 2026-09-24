# TSPCA

## Description

This repository includes the implementation of Tied Sparse PCA (TSPCA), a sparse PCA method in which the component weights are tied to the loadings (W = P) and a cardinality constraint fixes the number of nonzero weights. The method is implemented in two variants: a constrained approach, which imposes W = P exactly, and a penalized approach, which adds a penalty on the difference between W and P.

A simulation study compares TSPCA with sparse PCA methods that impose sparseness on either the weights (Zou et al.'s SPCA from the elasticnet package; GPower by Journée et al.) or the loadings (seafar). An illustration on the Big5 data shows TSPCA's advantages in explained variance and in stability of the zero/nonzero status of the weights. The repository contains the R functions that implement the method and the code for the simulation study and the illustration.

## Repository Structure

The repository is organised into two main folders:

-   Scripts/: Contains the R scripts for data generation and the functions required to run the methods.

    -   Psparse_Data.R : Script to obtain data with sparse loadings.

    -   Wsparse_Data.R : Script to obtain data with sparse weights.

    -   WPsparse_Data.R : Script to obtain data with sparse weights and loadings, including the symmetric condition W = P.

    -   TSPCA_Functions.R : Functions for the TSPCA method and for evaluating the results.

    -   GPower_Functions.R : Functions for the GPower method with a cardinality constraint.

-   Demo/ : Provides the scripts to run the methods on the simulated datasets and summarise the results.

    -   TSPCA_Sim.R : Runs TSPCA on all data-generation schemes and summarises the results per design cell, including the stability of the variable selection across replications.

    -   SPCA_Sim.R : The same for Zou et al.'s SPCA (elasticnet).

    -   GPower_Sim.R : The same for GPower.

    -   Seafar_Sim.R : The same for seafar.

    -   Illustration_Big5.qmd : Illustration on the Big5 data, comparing explained variance and resampling stability of TSPCA, SPCA, GPower and seafar.

## How to Run This Project

-   **Step 1:** Clone the repository

-   **Step 2:** Generate the data by running Scripts/Psparse_Data.R, Scripts/Wsparse_Data.R and Scripts/WPsparse_Data.R

-   **Step 3:** Move to the Demo folder and run the simulation scripts (TSPCA_Sim.R, SPCA_Sim.R, GPower_Sim.R, Seafar_Sim.R). Each script runs over all data-generation schemes in parallel.

## **Acknowledgements**

This project started as a Bachelor End Project, submitted in partial fulfillment of the requirements of the degree of Bachelor of Science at Eindhoven University of Technology and Tilburg University under the supervision of Prof. Dr. Katrijn Van Deun. This publication is part of the project SEM2.0 (with project numbers 406.22.GO.022 and VI.C.231.092) of the Open Competition and Talent research programs financed by the Dutch Research Council (NWO), awarded to Prof. Dr. Katrijn Van Deun.

## Authors

Bachelor Project by: Hetvi Chaniyara

Under the supervision of Prof. Dr. Katrijn Van Deun
