# FunCC
The FunCC package provides the funcc_biclust function that jointly performs clustering of rows and columns for a functional dataset. 
dakmapp

## Getting Started

The source code can be cloned or downloaded directly from github. An R studio project file is provided to open the project in RStudio.


## Installing

The package can be installed directly from github but devtools is required.

If devtools is not installed use the following comand to install it

install.packages('devtools') 
and then install the package

library(devtools)
install_github('mrtgalvani/FunCC')

Othewise the problem can be solved following the tutorial at: https://clang-omp.github.io/.

## Automatic tests

devtools::test()
## Example
data(funccdata)
res<-funcc_biclust(funccdata, ,delta=0.005,theta=1,number=500, alpha=0,beta=0,const_alpha=T,const_beta=T)
funcc_show_results(funccdata,res)

## Documentation

The R documentation can be found in the main directory in FunCC.pdf. 
