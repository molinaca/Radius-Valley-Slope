This repository contains all of the code and data files I used/compiled during my co-op term. It contains two main parts for each method that I used. The ABC folder will contain all the files that I used to
measure the occurence rates using Approximate Bayesian Computing. It includes scripts and notebooks for both the occurence rate method, and gaussian mixture method. The KDE folder contains a notebook used to create 
a kernel density estimate of planets from Berger et al. 2023. 

Packages needed:
- cosmoabc

The documentation for cosmoabc can be found at this link: https://pypi.org/project/cosmoabc/

Many of the functions in this were taken from Kunimoto et al. 2020, and can be found in the following repo: https://github.com/mkunimoto/Exo-Occurrence/blob/master/Codes/func_ExoOccurrence.py

**Contents:**

- `ABC/results_analysis.ipynb`: Notebook used to analyze data from both occurence rate and gaussian mixture method.
- `ABC/Occurrence_Rates/`
  - `my_make_input_file.py`: Creates an input_file to put into cosmoabc.
  - `job_utils.py`: Obtains information on period and radius bins given a job_id. `JOB_ID` must be defined in the environment as `JOB_ID=#`.
  - `radius_valley_func.py`: Contains all functions to create a forward model and distance function to input in cosmoabc.
  - `pick_file.py`: Run in directory with cosmoabc output to obtain final file numbers for each job_id.
- `ABC/Gaussian_Mixture/`
  - `make_input_gauss.py`: Create input_file to put into cosmoabc for the gaussian mixture method model.
  - `rvalley_gauss_func.py`: Contains all functions to create a forward model and distance function to input in cosmoabc.
- `KDE/kde_testing.ipynb`: Notebook that creates kde's for planets in the `GKTHCatalog_Table5.csv` table.

All data files used can be found in the Data folder. The `FGK_planets.dat` and `FGK_properties.dat` files were from the Kunimoto et al. 2020 Exo-Occurrence repo. The `GKTHCatalog_Table5.csv` table is from
Berger et al. 2023. 
      
**Key Notes**

- The occurence rate method is built to run multiple jobs (the codes are built around a single job_id). These are the codes that I ran on Cedar.
- I have still not figured out the issue with error from the Gaussian Mixture model.
- When running cosmoabc, I ran into some issues with where the packages were, this is the code to run it that worked for me: `python /home/molinaca/.local/bin/run_ABC.py -i input_file_${JOB_ID}.txt -f radius_valley_func.py`
