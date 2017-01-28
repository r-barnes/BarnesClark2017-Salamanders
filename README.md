BarnesClark2017-Salamanders
===========

**Title of Manuscript**:
65 Million Years of Change in Temperature and Topography Explain Evolutionary
History in Eastern North American Plethodontid Salamanders

**Authors**:
 * Richard Barnes (richard.barnes@berkeley.edu)
 * Adam Clark (adam.tclark@gmail.com)

**DOI Number of Manuscript**: TODO

**Code Repositories**
 * [Author's GitHub Repository](https://github.com/r-barnes/BarnesClark2017-Salamanders)



This repository contains source code for the salamander speciation model
described in the aforementioned manuscript.



Abstract
--------
For many taxa and systems, species richness peaks at mid-elevations. One
potential explanation for this pattern is that large-scale changes in climate
and geography have, over evolutionary time, selected for traits that are favored
under conditions found in contemporary mid-elevation regions. To test this
hypothesis, we used records of historical temperature and topographic changes
over the past 65 Myr to construct a general simulation model of Plethodontid
salamander evolution in eastern North America. We then explore possible
mechanisms constraining species to mid-elevation bands by using the model to
predict Plethodontid evolutionary history and contemporary geographic
distributions. Our results show that models which incorporate both temperature
and topographic changes are better able to predict these patterns, suggesting
that both processes may have played an important role in driving Plethodontid
evolution in the region. Additionally, our model (whose annotated source code is
included as a supplement) represents a proof of concept to encourage future work
that takes advantage of recent advances in computing power to combine models of
ecology, evolution, and earth history to better explain the abundance and
distribution of species over time.



Compilation
-----------

The program can be produced simply by running **make** in the root directory.
OpenMP is required for compilation.

Running `make` will produce an executable called `salamander.exe`.

Entering the `src` directory and running `make test` will generate `test.exe`,
which can be used to verify some parts of the code.



Running the Program
-------------------

`salamander.exe` can be run without arguments. In this mode it prints out the
format of its configuration file. An example configuration file is included at
`params/z_example_param.param`.

Running with a configuration file, as in

    ./salamander.exe params/z_example_param.param

generates several output files in the `output` directory.



Output Files
------------

The output files are in CSV format and contain self-documenting headers.
The `.tre` files contain the salamanders' phylogenies in Newick format.



Data Files
----------

The `data` directory contains the temperature time series referred to in the
manuscript. The file is ordered such that the most recent data points at at the
top. The intervals between lines are 1,000 years and the data series begins 65
Mya.
