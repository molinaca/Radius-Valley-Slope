#Cosmoabc will print out a data file for each iteration/particle system made. This file creates an array of the final iteration file numbers for each job right before it converged. 
#It will also print out the job numbers, this is useful for checking if all the jobs were successful.

import os
import re
from collections import defaultdict

# Define directory containing the files
directory = "/home/molinaca/Radius_Valley/Output/"  # Change this to your actual directory

# Regular expression to extract n and m from filenames
pattern = re.compile(r"rvalley_sims_(\d+)_(\d+)\.dat")

# Dictionary to store max m for each n
max_m_per_n = defaultdict(int)

# Loop through all files in the directory
for filename in os.listdir(directory):
    match = pattern.match(filename)
    if match:
        n, m = map(int, match.groups())  # Extract and convert to integers
        max_m_per_n[n] = max(max_m_per_n[n], m)  # Update max m for this n

# Convert to sorted dictionary for readability
max_m_per_n = dict(sorted(max_m_per_n.items()))

# Print results
#for n, max_m in max_m_per_n.items():
    #print(f"Bin {n}: Max m = {max_m}")

n_array = list(max_m_per_n.keys())
max_m_array = list(max_m_per_n.values())
print("n array:", n_array)
print("max m array:", max_m_array)