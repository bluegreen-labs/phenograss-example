# PhenoGrass Model Example

This is a working example of the PhenoGrass model as described in the Nature Climate Change paper: ["Productivity of North American grasslands is increased under future climate scenarios despite rising aridity"](http://www.nature.com/nclimate/journal/vaop/ncurrent/full/nclimate2942.html) by myself and co-authors.

### Figure 1

In particular this is a worked example which recreates Figure 1 of the journal article. It therefore includes a subset of the original data and uses the parameters and model framework as used throughout the rest of the publication.

Either download the zipped project or clone the project using:

  git clone https://github.com/khufkens/phenograss.git

Subsequently, run the Figure_1.r R script in the created directory:

On Linux you can call this script directly using:

  ./Figure_1.r

On OSX you can source the file in the R environment:

  source("Figure_1.r")

All model output and the generated PDF figure can be found in the output subdirectory within the root directory of the project. The result should look like the image below.

![](https://farm2.staticflickr.com/1524/26288199306_4b534c1202_o_d.png)

### Notes:

The code depends on a working gfortran compiler, which is freely available for all OS. However, you are on your own when it comes to rewriting the R code for Windows. The current code should run of the box on Linux and OSX
