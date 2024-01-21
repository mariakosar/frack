## About The Project

This project contains an extensive overlook of the data prepataion and modelling process for my dissertation chapter co-authored with Veronika Selezneva and Sergei Seleznev "Do financially constrained firms engage in opportunistic and risky behavior? Evidence from the oil and gas producing sector." We document opportunistic behavior of financially constrained firms in the oil and gas production sector. Our results show that firms exposed to debt shocks tend to increase the usage of experimental fracturing fluid combinations during well treatment jobs, a risky practice that can have significant impact on ultimate well recovery. The experimentation in response to the debt shock is delayed, which indicates that firms are not willing to take immediate risks. Firms with low number of debt covenants and firms with unsecured debt tend to experiment at a higher rate. 

The full (preliminary) article can be found here: https://www.overleaf.com/read/sfmmbhdqmghk#7b3695
 
### FracFocus data

FracFocus is a website of the Chemical Disclosure Registry of US. The site was created to provide the public access to reported chemicals used for hydraulic fracturing within the frack area.

1. Get a free FracFocus dataset at https://fracfocus.org/data-download.
   
2. Use `ffr.R` script to clean the FracFocus data:
    - check validity of CAS number. Citation: `is.cas_f` function is directly based on source code for `is.cas` function from `webchem` package: https://www.rdocumentation.org/packages/webchem/versions/1.1.2/topics/is.cas
    - replace invalid CAS with missing
    - remove leading zeros
    - sum similar CAS by UploadKey (single frack job), keep unique uploadkey and CAS
    - change dataset from long to wide format
    - create registry table
    - remove duplicated APINumber and JobStartDate
    - remove all entries with JobStartDate<31 May 2013 & >31 Dec 2020
    - merge registry with wide chemicals data
    - keep unique wells
    - exclude water, sand, and unknowns
    - clean data

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Jaccard index

The script jaccard.R computes Jaccard index for newly drilled oil and gas wells by comparing similarity of chemicals used for fracking well j to chemicals used for fracking well i. Pairwise similarities will be used for computing distances between well i and j in terms of proportions of chemicals used in fracking. Whenever well j is distant enough from a sufficient number of priorly drilled wells, we will call this well experimental well. 

We follow the methodology of [Fetter et al., 2018] to construct a Jaccard index as a measure of similarity. Pairwise similarity index between wells i and j is defined as:

![image](https://github.com/mariakosar/frack/assets/17361605/34427aa2-72bf-4f44-92ee-fd76cebc1564)

The main advantage of Jaccard index is that it is bounded by [0,1]. s_ij = 1 means that the fluids used in the fracturing process of wells i and j are identical, whereas s_ij = 0 implies that fluids used in wells i and j are completely different.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Selection of clustering algorithm

We need to identify whether chemical mixture of well i is dissimilar enough to wells that were developed prior to it. One way to do this is to use an unsupervised learning algorithm, such as density-based clustering algorithm (DBC). DBC identifies high-density regions (regions with large number of points within an ≤ e-ball with specified diameter) that are separated by regions of low density. The advantage of such
algorithm is that it can find clusters of arbitrary shape, and is capable of detecting outliers that do not fit any cluster ("noise points"), allowing to identify the experimental well ([Kriegel et al., 2011]).

#### DBSCAN

<p align="right">(<a href="#readme-top">back to top</a>)</p>

#### OPTICS
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>

#### LOF
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>

#### Heuristic clustering approach
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>

#### Processing of geo data in R
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>









