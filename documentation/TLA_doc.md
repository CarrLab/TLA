# Tumor Landscape Analysis (TLA)
 
Landscape ecology analysis methods for the investigation and characterization of digital histopathological data of tumor biopsies.

For more detailed explanations of these methods please refer to the publication: [Cisneros et al (2025)](http://www)


## Introduction

TLA implements different levels of spatial statistics methodologies in the study of tissues, perceived as cellular ecologies in the context of tumor development, and imaged by means of digital histopathologically. The data used for this analysis typically comes from a segmentation and classification process, which consist of image processing and machine learning algorithms that are capable of identifying tissue compartments or individual cells from a histopathological image, returning a raster map or a set of cell locations and cell type categories. In some cases, a combination format structure can include additional _ad hoc_ masks representing regional segmentations, also detected by image processing techniques, that identify specific tissue compartments of interest that can be used as a regional mask superimposed to the point-level data, like ducts or crypts, where neoplastic epithelial cells might form tumors or cellular niches. Similarly, the input data might consist only of regional segmentations indicating patches of tissue components, like tumor and stroma regions, that might be of interest to investigate in themselves.

## Point Statistics 

These are statistics calculated directly from the distribution of point events, representing cell locations and classifications, and specifically from the distributions of distances or aggregation behavior.

These statistics are done __BOTH__ globally (returning a single index for the whole study sample) and locally (returning local values evaluated in the vicinity of each location, and represented as raster profiles). 

In order to get local measures of spatial statistics, a landscape is typically segmented into quadrats for which each statistics is evaluated. In many cases, each quadrat is further segmented into sub-quadrats in order to perform the spatial calculations. The calculated index is then assigned to the whole quadrat. 
This TLA implementation use a convolution smoothing function:

$$y_{m,n} = \sum_i\sum_j x_{i,j}\cdot w_{m-i, n-j}$$

with $w$ the weighted kernel for the spatial convolution.
TLA uses a Gaussian kernel as a replacement for quadrats, and sub-quadrats, with bandwidths specified by the user. Thus Instead of segmenting the landscape into a grid to estimate spatial statistics, it estimates the statistical measure in each kernel around each pixel. 
This allows for calculating full resolution fields in a computationally efficient way.

### 1. First-order metrics:

These are conformation measurements that account for the intensity of events of each type. They involve __total numbers of cells of each type, total density of cells of each type, and fractions of the total of cell of each type__. In cases where an ad hoc filter mask is available, these values are accounted for each region of interest (as specified by said filter).

Local profiles of these measures are calculated by the use of a __Kernel Density Estimator (KDE)__, a spatial Gaussian convolution of discrete point events that estimates local, and smooth, point density values for each cell type, or combination of cell types.  

### 2. Colocalization index

Pairwise Morisita-Horn index between two classes is measured as the overlap of the spatial distributions:

$$M_{rt} = \frac{2 \cdot\sum_i{\left(n_i \cdot m_i\right)}}{\sum_i{(m_i)^2} + \sum_i{(n_i)^2}}$$

with the distributions of abundances of "ref" species _r_ and "test" _t_ species in each location _i_ of a spatial binning (i.e. histogram) given as:

$$
\begin{array}{lcr}
N_r & = & \sum_i{n_i} \\
N_t & = & \sum_i{m_i}
\end{array}
$$

__Note:__ because this calculation doesn't mix point abundances from different bin locations, the 2D spatial distributions can be represented as simple 1D histograms, or arrays, for the purpose of this calculation.  

On the other hand, this index can be calculated __locally__ in each region designated by the center _(x,y)_ as the overlap of __local__ spatial distributions:

$$M_{rt}(x,y) = \frac{2 \cdot\sum_i{\left(n_i(x,y) \cdot m_i(x,y)\right)}}{\sum_i{(m_i(x,y))^2} + \sum_i{(n_i(x,y))^2}}$$

with the distributions of abundances of "ref" species _r_ and "test" _t_ species in each sub-region _i_ of the _(x,y)_ region given as:

$$
\begin{array}{lcr}
N_r(x,y) & = & \sum_i{n_i(x,y)} \\
N_t(x,y) & = & \sum_i{m_i(x,y)}
\end{array}
$$

Globally or locally, this index measures the degree of overlap between the two corresponding spatial distributions, telling whether the two distributions are similar or not:

* _M_ ~ 0: ref and test cells are excluded from each other
* _M_ ~ 1: ref and test cells are similarly distributed in space

__Attention:__ this measure is very sensitive to low seeding.  Particularly when one or both distributions have empty bins. Namely, when the number of points of a particular cell type is significantly less than the number of bins (i.e. less than one cell per bin in average) it's statistically impossible to distinguish a discrete realizations of a non-uniform distribution from a low-seeded uniform distribution. Therefore for low-seeding cases, the $M$ would be overestimated and have a value with very low confidence.

This measure is symmetric, characterizing the prevalence and heterogeneity in colocalization of two cell types. All pair combinations are calculated.

### 3. Nearest Neighbor index

Measures the enrichment of nearest neighbors between pairs of cell type classes, "ref" cell type (blue dots) and the "test" cell type (red dots):
<div align="center">
<img src="nnind.png" width="150">
</div>

For each ref cell the distance to closest ref cell _d<sub>min</sub>(r)_ and the distance to the closest test cell _d<sub>min</sub>(t)_ are found and averaged, then the index is given by the ratio:

$$N_{rt} = \log\left(\frac{\langle d_{\text{min}}(t)\rangle}{\langle d_{\text{min}}(r)\rangle }\right)$$

In the case of the local profile, an identical calculation is performed in each individual region _(x,y)_ (which is also the region over which means are calculated): 

$$N_{rt}(x,y) = \log\left(\frac{\langle d_{\text{min}}(t)\rangle_{(x,y)}}{\langle d_{\text{min}}(r)\rangle_{(x,y)} } \right)$$

This measure has the properties:

* _N_ > 0: if ref and test cells are segregated from each other, meaning next to a ref cell is more likely to find another ref cell than a test cell.
* _N_ ~ 0: if ref and test cells are well mixed together, meaning the next to a ref cell is equally likely to find a ref cell or test cell
* _N_ < 0: if ref cells are individually infiltrated, meaning that next to a ref cell is more likely to find a test cell than a ref cell. This will happen if test cell are mixed in but much less abundant.

This measure is different from the colocalization score in that it captures the character of cell-cell closeness between different cell types, typically imposing a weaker condition than overall spatial distribution similarity. Namely, if two cell types are "colocalized" in the "Morisita-Horn index sense" (_M_ ~ 1) then typically _N_ ~ 0, but the converse is not necessarily true: two cell types can be "proximal to each other", _i.e._  _N_ ~ 0, and still have very different spatial distributions, specially if the abundances of the two types is very different. Therefore the spatial profiles of this factor characterize the heterogeneity in mixing of cells in terms of their typical observed proximity.

This is an asymmetric pairwise measure that is sensitive to relative abundance of cells.

### 4. Ripley’s H Index 

This factor measures relative clustering of points as a function of scale of the scope. In this case, for a given radius $d$ and cell types "ref" (blue dots) and "test" (red dots):
<div align="center">
<img src="rhind.png" width="150">
</div>

For each ref cell the number of test cells $I_{rt}(d)$ inside a radius $d$ from it is calculated. The mean across all ref cells normalized by the density $\lambda_t$ of test cells gives the Ripley's $K$ function:

$$ K_{rt} = \frac{1}{\lambda_{\text{t}}}\langle I_{\text{rt}}(d)\rangle$$

If the distribution of test points is homogeneous, i.e. they satisfy a null hypothesis of __Complete Spatial Randomness (CSR)__, the expected value of $I$ should approach $E = A\lambda_t$ with $A = \pi d^2$ the area of the circle. Therefore $K$ is given by the proportion of observed to expected points in the circle. The more practical $H$ function is defined as:

$$H_{rt} = \log\left(\sqrt{\frac{K_{rt}}{\pi d^2}}\right)$$

For the local profile of this factor, the calculations are done in each region _(x,y)_ as:

$$ K_{rt}(x,y) = \frac{1}{\lambda_{\text{t}}}\langle I_{\text{rt}}(d)\rangle_{(x,y)}$$
$$H_{rt}(x,y) = \log\left(\sqrt{\frac{K_{rt}(x,y)}{\pi d^2}}\right)$$

This is a measure of the level of clustering of test cells around ref cells at the scale _d_ and relative to a CSR null hypothesis. Typically, this measure is used to assess the structure of a spatial distribution, generating a curve as a function of _d_ that characterize how the clumping changes with the distance, leading to estimation of natural scales and correlation distances in the system. Given that TLA generates local estimates of each factor, it uses only one arbitrary value, _d_ equals sub-region size, to evaluate the Ripley's H function as a spatial factor. It has the following properties:

* _H_ > 0: if test cells cluster around ref cells
* _H_ ~ 0: if ref and test cells are mixed uniformly
* _H_ < 0: if test cells are excluded from ref cells

This is yet another factor that describes the mixing of two class types in a slightly different way. This measure refers to the level of clumping at a specific scale _d_. It is a way related to the Nearest Neighbor index, more stringent in the continuous aspect of the cell's spatial distribution but yet limited to a length scale. Therefore, this is an intermediate approach in terms of spatial distribution comparisons. 

This is not only an asymmetric measure, but also not an identity. Applying the measure on ref cells with themselves (i.e. test cells are the same as ref cells) does not necessarily return a trivial identity value, because the comparison is against a CSR hypothesis rather than the observed distribution. The identity figure is then a measure of the degree of clumping of ref cells around themselves, thus producing a measure of effective clustering. Pairwise comparisons give a measure of relative clumping, or clustering, of test cells __around__ ref cells. 


### 5. Getis-Ord Gi* score and HOT index

This is a standard measure of enrichment of cell abundance with respect to the total study area (landscape). It consist of a general inferential _Z_ statistic, tested in the context of the null hypothesis. See [ArcGIS documentation](https://pro.arcgis.com/en/pro-app/2.8/tool-reference/spatial-statistics/h-how-high-low-clustering-getis-ord-general-g-spat.htm) for mathematical expressions and definitions. Thus, this is by construction a local profile metric that has no global form. 

The null Complete Spatial Randomness (CSR) hypothesis states that there is no spatial clustering of feature values. When the corresponding p-value is statistically significant the null hypothesis can be rejected in one of two ways according to the sign of the z-score: if positive, the observed index is larger than expected, indicating that high values for the attribute (cell abundance) are clustered in the study area. If the z-score value is negative, the observed index is smaller than expected indicating that low values for the attribute (cell abundance, hence depletion) are clustered in the study area.

Along with a spatial factor for the z-scores, TLA generate a HOT factor consisting of:

* HOT ~ 1: if region _(x,y)_ has more cells than expected
* HOT ~ 0: if region _(x,y)_ has a expected density (null hypothesis cannot be rejected)
* HOT ~ -1: if region _(x,y)_ has less cells than expected

These spatial profiles characterize hotspots in relation to the __whole study area__ (sample landscape), indicating regions of abundance or depletion for each cell type. Overlaps for different classes can be generated to detects mutual hot/cold spotting in different combinations. 


## Landscape Ecology

Landscape ecology is the study of interactions and relationships between living organisms and the environment they inhabit. Such environment is defined as a __landscape__ (spread in space and time) occupied by different species of organisms, and its description entails details on spatial distributions, cohabitation, population dynamics, mobility and other ecological questions.

A __landscape mosaic__ consist of spatial locations and categorical classification of all ecological entities of interest. Typically such data is recorded as a __categorical raster image__ (a 2D array with discrete value labels representing surveyed points in a particular moment in time and location, with different categories representing species), or an equivalent long-form table of individuals with their corresponding coordinates and category values. 

Landscape metrics algorithms have been implemented in packages such as [FRAGSTATS](http://www.umass.edu/landeco/research/fragstats/fragstats.html), [landscapemetrics](https://r-spatialecology.github.io/landscapemetrics/) (R) or [pylandstats](https://pylandstats.readthedocs.io/en/latest/index.html) (python) support such raster spatial objects. Because these algorithms work with categorical data, each cell in the array is assigned to a discrete class identifier. Therefore there is an inherent spatial resolution of the data given by the discrete location values. We will refer to this resolution as the "pixel resolution" of the data.

__Landscape metrics__ are measurements that quantify physical characteristics of landscape mosaics that relate to ecological processes. These tools help characterize a landscape, primary describing its composition and spatial configuration: 

- The __composition__ of a landscape accounts for how much of the landscape, or a specified region of it, is covered by a certain category type or any measure of density, frequency or intensity of effects.
- The __configuration__ describes the spatial distribution or complex arrangements of the different categories as an effect of the interactions between conforming elements of the landscape. 

Additionally, landscape metrics can be calculated for three different levels of information scopes:

1.	__Patch level metrics:__ a patch is defined as a region of neighboring locations belonging to the same class, typically using Moore's or von Neumann's neighborhood rules. Patch level metrics are calculated for each patch in the landscape.
2. __Class level metrics__: returns a summary value for all patches aggregated by type class. The output is typically some statistics of patch-level metrics across all patches in each class (e.g. a sum or mean). 
3. __Landscape level metrics__: returns a single value describing a global property of the landscape. This is typically a statistics of metrics of lower levels aggregated by patches and/or classes or a multivariate metric calculated globally (like a diversity or entropy measure). 

There are different classes of landscape metrics:

1.	__Area and edge metrics__: describe the aggregated area and length of the edge of patches or classes. The __edge__ is defined as the border perimeter of patches. These metrics mainly characterize the composition of the landscape and evaluate dominance or rareness of classes.
2.	__Shape metrics__: describe the shape of patches, typically in terms of the relationship between their area and perimeter but also including or other metrics that describe the overall geometry of each patch (like fractal dimension, roundness, etc). 
3.	__Core area metrics__: describe the area of the fraction of each patch that is not an edge, providing information about areas that are not influenced by neighboring patches of a different class.
4.	__Contrast metrics__: describe the magnitude of the difference between adjacent patch types with respect to one or more ecological attributes. 
5.	__Aggregation metrics__: describe the level of clumpiness of patches of the same class, providing information on whether patches or a certain class tend to be aggregated in space or isolated. These metrics describe the spatial configuration of the landscape.
6.	__Diversity metrics__: available on the landscape level, these metrics describe the abundance and dominance/rareness of classes and show the diversity of classes in the landscape.
7.	__Complexity metrics__: provide information theory-based measures, like entropy and mutual information, to characterize patches of given classes.  

### Metric statistics

Toolboxes like __pylandstats__ (which is implemented within TLA) feature six specific distribution metrics for each patch-level metric, consisting of statistical aggregation of the values computed for each patch of a class or the whole landscape.

In what follows we refer to the following notation:  

- __a<sub>i,j</sub>, p<sub>i,j</sub>, h<sub>i,j</sub>__ represent the area (__a__), perimeter (__p__), and distance to the nearest neighboring patch of the same class (__h__) of the patch __j__ of class __i__.
- __e<sub>i,k</sub>, g<sub>i,k</sub>__ represent the total edge (__e__) and number of pixel adjacencies (__g__) between classes __i__ and __k__.
- __A, N, E__ represent the totals of area (__A__), the number of patches (__N__) and edge (__E__) of the landscape

Six basic distribution metrics are calculated across patches (and done at the class-level or landscape-level):

1.	__Mean__: specified by the suffix `_mn` to the method name, e.g. `area_mn`. 
2.	__Area-weighted mean__, specified by the suffix `_am` to the method name, e.g. `area_am`. This is the mean value weighted by the patch size.
3.	__Median__, specified by the suffix `_md` to the method name, , e.g. `area_md`. 
4.	__Range__, specified by the suffix `_ra` to the method name,  e.g. `area_ra`.
5.	__Standard deviation__, specified by the suffix `_sd` to the method name,  e.g. `area_sd`
6.	__Coefficient of variation__, (or variance) specified by the suffix `_cv` to the method name,  e.g. `area_cv`


## Spatial Stratified Heterogeneity

Spatial heterogeneity refers to non-uniform distribution of spatial factors within an area. Environments can have a variety of habitats, such as different topographies and ecological diversities capable of accommodating more or less interacting species: when organisms can finely partition a landscape into distinct suitable habitats (or niches), more species can coexist with minimal competition.

Spatial stratified heterogeneity refers to a situation in which landscape strata, such as ecological zones or classes, present distinct profiles for spatial factors such that within strata variance is less than the between strata variance. This indicates that the stratification is somewhat predictive or well correlated to the intensity profile of the spatial factor.

Please refer to documentation for the R package [geodetector](https://cran.r-project.org/web/packages/geodetector/vignettes/geodetector.html) for more details, definitions and examples.

### 1. SSH Factor detector

The factor detector q-statistic measures the Spatial Stratified Heterogeneity (SSH) of a spatial factor _Y_ in relation to a categorical variable _X_ (strata). This is also known as the determinant power of a covariate _X_ of _Y_. Outputs include the following statistics:

* __q-statistic__: The q-statistic measures the degree of SSH: 
	- __q~1__ indicates a strong stratification: small within-strata
variance and/or large between-strata variance. Thus a strong association between the explanatory variable and the explained variable (ie. strata categories explain the data)
	- __q~0__ no stratification: within-strata variance is large and/or between-strata variance is small. Thus there is no relationship between the strata categories and the data.
* __F-statistic__: F-value, assuming a random non-central F-distribution
* __p-value__: Prob that the q-value is observed by random chance. The null hypothesis is defined as absence of within-stratum heterogeneity (__q~0__):
	- __H<sub>0</sub>__: there is no SSH (stratification is not significant), thus within and between strata heterogeneity are similar.
	- __H<sub>1</sub>__: there is SSH (stratification is significant), thus within-strata heterogeneity is significantly smaller than between-strata heterogeneity.

### 2. SSH Interaction detector

The interaction detector function reveals whether risk variables
_{X<sub>1</sub>, X<sub>2</sub>}_ have an interactive influence on a factor _Y_. Outputs a table for the interactive q-statistics between variables, accounted by measuring the effect of merging _X<sub>1</sub>_ and _X<sub>2</sub>_ together in a new category defined by their union, giving the possible outcomes:

* __"equivalent"__ if q(X1∩X2) = q(X1) = q(X2)
* __"weaken"__ if q(X1∩X2) < q(X1) + q(X2)
* __"weaken, nonlinear"__ if q(X1∩X2) < q(X1) and q(X2)   
* __"max(q(X1),q(x2)) weaken (uni-)"__ if q(X1∩X2) < max(q(X1),q(x2))
* __"max(q(X1),q(x2)) weaken; min(q(X1),q(x2)) enhance"__ if min(q(X1),q(x2)) < q(X1∩X2) < max(q(X1),q(x2))
* __"min(q(X1),q(x2)) enhance (uni-)"__ if min(q(X1),q(x2)) < q(X1∩X2) 
* __"independent"__ if q(X1∩X2) = q(X1) + q(X2)
* __"enhance, nonlinear"__ if q(X1∩X2) > q(X1) + q(X2)
* __"enhance, bi-"__ if q(X1∩X2) > q(X1) and q(X2)
   
### 3. SSH Risk detector

This function calculates the average values in each stratum of the explanatory variable _X_, and reports if a significant difference between any two strata levels exists, indicating that the factor is a risk factor for the stratification structure. It outputs means of explained variable in each stratum and the t-test for differences between every pair of strata (with the corresponding multiple comparison correction for p-values)

### 4. SSH Ecological detector

This function identifies the impact of differences between two risk factors _X<sub>1</sub>_ and _X<sub>2</sub>_  returning the significance test of impact difference between them.

The probability value is calculated from the positive tail of the cumulative F-statistic given by the ratio of the two factors individual F-statistics. The significance of this measure indicates that the two stratifications _X<sub>1</sub>_ and _X<sub>2</sub>_ are statistically distinct in terms of risk.

## References:

1. Mcgarigal, K., Cushman, S., & Ene, E. (2012). FRAGSTATS v4: Spatial Pattern Analysis Program for Categorical and Continuous Maps. Retrieved from http://www.umass.edu/landeco/research/fragstats/fragstats.html
2. Hesselbarth, M. H. K., Sciaini, M., With, K. A., Wiegand, K., & Nowosad, J. (2019). landscapemetrics: an open-source R tool to calculate landscape metrics. Ecography, 42(10), 1648–1657. https://doi.org/10.1111/ecog.04617
3. Nowosad, J., & Stepinski, T. F. (2019). Information theory as a consistent framework for quantification and classification of landscape patterns. Landscape Ecology, 34(9), 2091–2101. https://doi.org/10.1007/s10980-019-00830-x
4. Wolda H. Similarity indices, sample size and diversity. Oecologia. 1981 Sep;50(3):296-302.
5. Magurran A.E. (2005) Biological diversity. Curr Biol 15:R116-8
6. Rempala G.A., Seweryn M. (2013) Methods for diversity and overlap analysis in T-cell receptor populations. J Math Biol 67:1339-68
4. Altman, Naomi S. (1992). "An introduction to kernel and nearest-neighbor nonparametric regression". The American Statistician. 46 (3): 175–185.
4. Everitt, Brian S.; Landau, Sabine; Leese, Morven; and Stahl, Daniel (2011) "Miscellaneous Clustering Methods", in Cluster Analysis, 5th Edition, John Wiley & Sons, Ltd., Chichester, UK
4. Bosch, M. (2019). PyLandStats: An open-source Pythonic library to compute landscape metrics. BioRxiv, (October), 715052. https://doi.org/10.1101/715052
5. Ripley, B.D. (1976). "The second-order analysis of stationary point processes". Journal of Applied Probability. 13 (2): 255–266. doi:10.2307/3212829. JSTOR 3212829.
6. Dixon, Philip M. (2002). "Ripley's K function". In El-Shaarawi, Abdel H.; Piegorsch, Walter W. (eds.). Encyclopedia of Environmetrics. John Wiley & Sons. pp. 1796–1803. ISBN 978-0-471-89997-6.
5. Getis, Arthur, and J. K. Ord. "The Analysis of Spatial Association by Use of Distance Statistics." Geographical Analysis 24, no. 3. 1992.
6. Mitchell, Andy. The ESRI Guide to GIS Analysis, Volume 2. ESRI Press, 2005.
7. Wang JF, Li XH, Christakos G, Liao YL, Zhang T, Gu X, Zheng XY.
Geographical detectors-based health risk assessment and its application in the neural tube defects study of the Heshun Region, China. International Journal of Geographical. Information Science, 2010, 24(1): 107-127.
7. Wang JF, Zhang TL, Fu BJ. 2016. A measure of spatial stratified heterogeneity. Ecological Indicators 67: 250-256.
8. Wang JF, Xu CD. Geodetector:Principle and prospective. Geographica
Sinica, 2017, 72(1):116-134.
8. Jiang B. 2015. Geospatial analysis requires a different way of thinking: The problem of spatial heterogeneity. GeoJournal 80(1), 1-13.
9. Song, Y and Wu, P (2021). “An interactive detector for spatial associations”. International Journal of
Geographical Information Science. doi:10.1080/13658816.2021.1882680

---


