# Data dictionary

Date: the timestamp of the data

Load: the electricity loads

Forecast_RTE: forecasted loads from Réseau de Transport d'Electricité 

WeekDays: the day of the week

BH: public holidy

tod: index for 48 half-hours from 0-47

Year: year, from 2012-2019

Month: month, from 1-12

Load.48: the lag 48 load, i.e. the load at the same hour at the previous day

Load.336: the lag 336 load, i.e. the load at the same hour at the previous 7-th day

Temp: temporature 

Temp_s95: smoothed temperature with smooth rate 0.95

Temp_s99: smoothed temperature with smooth rate 0.99

Temp_s95_min: minimum smoothed temperature with smooth rate 0.95

Temp_s95_max: maximum smoothed temperature with smooth rate 0.95

Temp_s99_min: minimum smoothed temperature with smooth rate 0.99

Temp_s99_max: maximum smoothed temperature with smooth rate 0.99

Summer_break: 10 means is the summer break period, 0 means otherwise

Christmas_break: 20 means is the Christmas break period, 0 means otherwise

WeekEnd: TRUE means the day is in weekends, FALSE means otherwise
