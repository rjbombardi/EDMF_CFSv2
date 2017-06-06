# EDMF_CFSv2
The Eddy-Diffusivity Mass-Flux (EDMF) Boundary Layer Parameterization ready for implementation into the CFSv2

The bcltrigger.f is a subroutine that calculates the Heathed Condensation Framework (HCF; Tawfik et al. 2014,2015a,b) trigger as published in Bombardi et al. (2016).
The bcltrigger module also contains the namelist variables used to call the EDMF parameterization.

The files compns_v.f and the gbphys_v.f are modified versions of the original subroutines used in the Climate Forecast System (CFSv2) model adaped to work

gbphys_v is the subroutine that calls the EDMF parameterization.



Bombardi RJ, Tawfik AB, Manganello JV, Marx L, Shin C-S, Halder S, Schneider EK, Dirmeyer PA, Kinter JL (2016) The heated condensation framework as a convective trigger in the NCEP climate forecast system version 2. JAMES, DOI: 10.1002/2016MS000668

Tawfik, A. B., P. A. Dirmeyer, and J. A. Santanello,  2015a: The heated condensation framework. Part I: Description and Southern Great Plains case
study, J. Hydrometeorol., 16, 1929-1945, DOI 10.1175/JHM-D-14-0117.1.

Tawfik, A. B., P. A. Dirmeyer, and J. A. Santanello, 2015b: The heated condensation framework. Part II: Climatological behavior of convective
initiation and land-atmosphere coupling over the Conterminous United States, J. Hydrometeorol., 16, 1946-1961, DOI:10.1175/JHM-D-14-0118.1.

