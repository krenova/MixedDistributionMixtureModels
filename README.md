# Mixed Distribution Mixture Models (MDMM)
This is a mixed distributed clustering model that uses the Gaussian Distribution and Multinomial Distribution to segment mixed typed data. For example, datapoints with a mixture of numerical, categorical or multinomial. 

Advantage of using mixed distributions to cluster mixed typed data is better segmentation due to the use of appropriate distributions to capture the respective data type characteristics. On a more fundamental level, such a clustering model also circumvent problems with using distance measure algorithms, for example [K-Means](https://en.wikipedia.org/wiki/K-means_clustering), to cluster categorical and multinomial data. While there are excellent algorithms out there that are designed to handle categorical and multinomial data, see [ROCK](https://www.google.com.sg/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwin9uHs98XSAhVJto8KHQAgCdwQFggYMAA&url=https%3A%2F%2Fwww.cis.upenn.edu%2F~sudipto%2Fmypapers%2Fcategorical.pdf&usg=AFQjCNGvLw8baJgnLQBMUQ3kfUckNuct4A&bvm=bv.149093890,d.c2I) clustering or [other variants of mixture models](https://en.wikipedia.org/wiki/Mixture_model), these algorithms do not work well with numerical data. MDMM is therefore designed to address these common problems faced when working with mixed data type. As MDMM is written using R's C++ api, [**Rcpp**](https://www.google.com.sg/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&cad=rja&uact=8&ved=0ahUKEwjhoPDb-cXSAhVGq48KHTpOAQUQFgghMAI&url=https%3A%2F%2Fcran.r-project.org%2Fpackage%3DRcpp&usg=AFQjCNGgJ-3UpTLaAHqVXWPK5RMdFSVKHw&bvm=bv.149093890,d.c2I), and uses the [Expectation-Maximization](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm) (EM) algorithm to infer model parameters, the code runs relatively fast and should have no problems handling data of a million data points.

## Codes
The main codes are:

1. [**mdmmCpp.R**](https://github.com/krenova/MixedDistributionMixtureModels/blob/master/lib/mdmmCpp.R)
2. [**mdmmCore.cpp**](https://github.com/krenova/MixedDistributionMixtureModels/blob/master/lib/mdmmCore.cpp)
	
where mdmmCpp.R is the file that your R script would source from and _**mdmmCore.cpp**_ is the C++ code that _**mdmmCpp.R**_ is built on.

For a demonstration on the use of the clustering function and also compare the speeds for an equivalent code written in R, do refer to the following jupyter notebook:

* [**demo.ipynb**](https://github.com/krenova/MixedDistributionMixtureModels/blob/master/demo.ipynb)


## Issues
As would be expected of an EM algorithm, the log-likelihood should be monotonically increasing. However, there are instances where the log-likelihood dips which goes against EM theory. Work is needed to determine if the cause of the dips is due to numerical overflows or bugs.



## Reading References
[1] http://datascience.stackexchange.com/questions/22/k-means-clustering-for-mixed-numeric-and-categorical-data