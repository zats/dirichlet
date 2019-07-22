# A new statistical method to analyze Morris Water Maze data using Dirichlet distribution

This is a python package forked from Eric Sue's [dirichlet](https://github.com/ericsuh/dirichlet) package, adding a test of uniformity for a single data set (ie equal fractions in each category) and plotting capability.

This was intended for the particular case of Morris Water Maze data as presented in Maugard and Doux, 2018. The test includes an approximate Bartlett-type correction for small samples as the test is based on a likelihood-ratio, which is only accurate for large samples.

## Use 

Download the repository as a zip file and unzip it.

### With python

 We provide a jupyter notebook in the `example` subfolder to perform the test and produce plots, see[example.ipynb](https://github.com/xuod/dirichlet/blob/master/example/example.ipynb) (if you have never used jupyter, see the [jupyter website](http://jupyter.org/)).

<img width="20%" height="20%" src="https://github.com/xuod/dirichlet/blob/master/example/3Tg.png"> <img width="20%" height="20%" src="https://github.com/xuod/dirichlet/blob/master/example/wt.png">

### With R

You can import the python `dirichlet` package using the `reticulate` package for R. Here's an example from the R console:
```
library(reticulate)
setwd('/path/to/dirichlet/example')
data_3tg<-data.matrix(read.csv('3Tg.csv'))
dirichlet <- import_from_path("dirichlet", path="../")
dirichlet$test_uniform(data_3tg, label='3Tg', do_MWM_correction=TRUE, verbose=TRUE)
```
The data should be a 2D matrix with samples in rows. The `plot` function can also be used but the plot is not as good.

