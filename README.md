A new statistical method to analyze Morris Water Maze data using Dirichlet distribution
=======================================================================================

This is a python package forked from Eric Sue's [dirichlet](https://github.com/ericsuh/dirichlet) package, adding a test of uniformity for a single data set (ie equal fractions in each category) and plotting capability. This was intended for the particular case of Morris Water Maze data as presented in Maugard and Doux, 2018. The test includes an approximate Bartlett-type correction for small samples as the test is based on a likelihood-ratio, which is only accurate for large samples.

Use 
---

Download the repository as a zip file and unzip it. Then see the [example](https://github.com/xuod/dirichlet/blob/master/example/example.ipynb) jupyter notebook to run the test and produce plots.

<img width="20%" height="20%" src="https://github.com/xuod/dirichlet/blob/master/example/3Tg.png"> <img width="20%" height="20%" src="https://github.com/xuod/dirichlet/blob/master/example/wt.png">
