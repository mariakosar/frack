## About The Project

This project contains an extensive overlook of the data prepataion and modelling process for my dissertation chapter co-authored with Veronika Selezneva and Sergei Seleznev "Do financially constrained firms engage in opportunistic and risky behavior? Evidence from the oil and gas producing sector." We document opportunistic behavior of financially constrained firms in the oil and gas production sector. Our results show that firms exposed to debt shocks tend to increase the usage of experimental fracturing fluid combinations during well treatment jobs, a risky practice that can have significant impact on ultimate well recovery. The experimentation in response to the debt shock is delayed, which indicates that firms are not willing to take immediate risks. Firms with low number of debt covenants and firms with unsecured debt tend to experiment at a higher rate. The full (preliminary) article can be found here: https://www.overleaf.com/read/sfmmbhdqmghk#7b3695
 
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
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>

### DBSCAN
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>

### OPTICS
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>

### LOF
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Heuristic clustering approach
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Processing of geo data in R
TBA
<p align="right">(<a href="#readme-top">back to top</a>)</p>









