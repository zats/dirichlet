A new statistical method to analyze Morris Water Maze data using Dirichlet distribution
=======================================================================================

This is a python package forked from Eric Sue's [dirichlet](https://github.com/ericsuh/dirichlet) package, adding a test of uniformity for a single data set (ie equal fractions in each category) and plotting capability. This was intended for the particular case of Morris Water Maze data as presented in Maugard and Doux, 2018. The test includes an approximate Bartlett-type correction for small samples as the test is based on a likelihood-ratio, which is only accurate for large samples.

Installation
------------

    pip install git+https://github.com/xuod/dirichlet.git

Use 
---

See the [example](https://github.com/xuod/dirichlet/blob/master/example/example.ipynb) jupyter notebook to run the test and produce plots.

![plop](https://github.com/xuod/dirichlet/blob/master/example/3Tg.png "3Tg.png")
