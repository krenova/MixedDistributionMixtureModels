# MixedDistributionMixtureModels
This is a mixed distributed clustering model that uses the Gaussian Distribution and Multinomial Distribution to segment mixed typed data. For example, datapoints with a mixture of numerical, categorical or multinomial.
Advantage of using mixed distributions to cluster mixed typed data is better segmentation due to the use of appropriate distributions to capture the respective data type characteristics. 

## Main codes
The main codes are

    1) mdmmCpp.R
    2) mdmmCore.cpp
	
where mdmmCpp.R is the file that your R script would source from and mdmmCore.cpp is the C++ code that mdmmCpp.R is built on.

