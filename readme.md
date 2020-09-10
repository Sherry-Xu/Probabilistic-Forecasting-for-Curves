# Probabilistic Forecasting for Daily Electricity Loads and Quantiles for Curve-to-Curve Regression

## Data 
### Abstract
The real data set consists of French daily electricity load data and temperature data from January 1, 2012 to December 31, 2019. The French electricity consumption data are collected from the website of the system operator RTE (Réseau de Transport d'Electricité): https://opendata.rte-france.com) at a temporal resolution of every half-hour (i.e. 48 points on each day). We obtained data from 96 meteostations in France from the website of the French weather forecaster Météo-France(https://donneespubliques.meteofrance.fr/). Temperature data are provided at a three hours resolution and interpolated with natural cubic splines at a half-hour resolution.

### Availability
The real data used in the manuscript can be find in the data directory with a description of the data dictionary.


## Code
The code directory contians all the necessary files to conduct simulation experiments and real data analysis.

* PreditInt_F.R contains all the functions that are used in SimulationExperiments.R and RealDataAnalysis.R

* SimulationExperiments.R is for simulation experiments using the settings specified in parameters.R

* RealDataAnalysis.R is for probabilistic forecasting for electricity load curves

* Alternative_GAM_ARMAerrorcorrection_Recursive.R, Alternative_SAR_ARMAerrorcorrection_Recursive.R and Alternative_SARX_ARMAerrorcorrection_Recursive.R are the codes for the three alternative models used in the paper


## Reproducibility workflow
* For simulations, first run the SimulationExperiments.R by specify the setting in parameter.R. After obtaining the results, use SimulationResultsAnalysis.ipynb to produce the figures and tables.

* For Real data anlaysis, the results of the proposed model can be obtained by runing RealDataAnalysis.R. Results from the three alternative models can be obtained from the other three code scripts for alternative models.


## Results illustration
An illustration of the out of sample forecasting for the year 2019 can be find in https://www.dropbox.com/s/hi274jlu8bx4tnl/OutofSampleForecasts2019.mp4?dl=0 
