# Drosophila larvae -
## A package for analyzing spiking/bursting in drosophila larvae motor neurons.

This package has grown organically and most of the time it lives in Python. Spike detection is still done in Matlab. Analysis has flowed like this:

* Get all the .abf files of the desired type ('GapFree', '10^-6 Pilo', etc); this can be done with getExtFromDF in the main package, _drophila.py_
* Pass that to _batchWriteSpikes.m_, which uses a moving window to retrieve all spikes from every file listed
* Every .abf file now has a .csv file that includes spike properties for every spike time
* If burst detection is desired, the main package can also find bursts using inter-spike intervals and Markov chain Monte Carlo simulation -- this is not perfect, but gets about 90% of the spikes. Burst inclusion (1=burst, 0=tonic) is appended to the trace's .abf file
* Specific burst features can be extracted if desired using tools in the main function including burst times, inter-burst intervals, and coefficients of variation for each of these

## Burst detection

Some details for the burst detection algorithm can be found in the document *burst_tonic.pdf*  in this repo.
