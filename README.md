# window_adjust.py: a program to adjust genomic window sizes along with average 
# window values.


## Overview 

**Goal**: My need for this script came about because the program ReLERNN produces a file of recombination
rates for genomic windows of variable lengths. For comparison with other genome stats, I wanted to adjust the 
recombination windows to the same length as the windows for all my other stats (in my case, 10 kb windows).

Thus, this program checks for overlap between the user-specified windows and the ReLERNN windows to appropriately
calculate a weighted average of recombination rate for the new 10 kb windows.

**Input**: 
- ReLERNN output file with five columns: chrom, start, end, nSites, and recombination rate
- genomic windows file (.bed)

**Output**:
- One file (*recomb_rate_by_10kb_windows.csv*) with the newly-calculated recombination rate averages for 10 kb windows.

## Command to run example file

`python window_adjust.py --dataFile xiph.hardfilt.9indTapajos.hiqual.RMexclude.complete.PREDICT.txt --windowFile windows.bed`

