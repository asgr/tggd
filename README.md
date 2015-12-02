tggd
====
This small package provides the standard R distribution functions for the Truncated Generalised Gamma Distribution  (see Murray, Robotham, Power (2016) [MRP]).

Included are the density, CDF, quantile function and random generator, both for the real-space TGGD and log-space TGGD. 

Installation
------------
The simplest way to install the package is to use ``devtools``. 

First, install the only dependency, ``gsl``: ``install.packages("gsl")``.

Next, install ``devtools``: ``install.packages("devtools")`` (note there are further dependencies of this package).

Finally, install this package: ``install_github('asgr/tggd')``.

Usage
-----
Usage details are found by accessing the help files in the package, but are identical to the standard R distribution functions (eg. ``rnorm`` or `dpoiss``). Examples are also found in these files.

To understand *why* and *when* to use these functions, please refer to MRP. 

We note in particular that the log-space forms of the functions are defined as the distributions of the log of the variates from the real-space distribution with identical parameters. This is the correct form to use if performing Bayesian analyses in which the values of the variates are uncertain with priors defined in log-space (eg. lognormal distributions).

Authors
-------
* Aaron Robotham
* Steven Murray

Version
-------
0.1
