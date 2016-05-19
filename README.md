spenra
======

What is spenra?
---------------

<b>spenra</b> is a package for computing the specific differential entropy rate of a continuous-valued time series. It falls into a similar class of statistics as Approximate Entropy (ApEn), Sample Entropy (SampEn), and Permutation Entropy, and provides an estimate of the intrinsic uncertainty associated with a time series. **spenra** uses modern non-parametrics to tune itself to the geometry and memory of the time series.

The accompanying paper is [here](http://www.mdpi.com/1099-4300/18/5/190). It provides several worked examples.

Installation
------------

The package can be installed with

	devtools::install_github("ddarmon/spenra")

After installation, the package can be loaded into R using

	library(spenra)

Using spenra
------------

See the files <tt>demo-condstat.R</tt> and <tt>demo-lorenz.R</tt> for example usage of **spenra**.

If you use **spenra** in a scientific publication, please cite

Darmon, D.	Specific Differential Entropy Rate Estimation for Continuous-Valued Time Series. *Entropy* **2016**, *18*, 190.