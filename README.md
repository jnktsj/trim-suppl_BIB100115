# Read quality trimming scripts

Read quality triming methods can be grouped into two main classes,
running-sum and window-based. To test the performance of read trimming
methods, we implemented the following trimming algorithms:

### `mott-trim.py`: Mott trimming (running-sum)
The method starts from the 3´-end of each read, subtracts a preset
cutoff quality score from the quality score at each position and adds
the remainder to a cumulative score at the position. The 3´ portion of
the read starting from the position with the minimum cumulative score
is trimmed.

### `dynamic-trim.py`: Dynamic trimming (window-based)
The method searches for the longest stretch of positions (window) in
each read such that the quality scores of each position in the window
exceed a preset threshold.

### `simple-trim.py`: **Simple trimming**
In addition to testing the two main trimming algorithm, we also
teseted simple trimming.  The method scans from the 5´ end of each
read. As soon as it detects a position with quality scores below a
preset threshold, it discards this position and the remaining
positions at the 3´-end of the read.

The trimming scripts are used in "Evaluation of pre-processing,
mapping, and post-processing algorithms for analyzing whole genome
bisulfite sequencing data".