# COVID19_vaccination-MPM
Modelling of regional vaccination strategies using a metapopulation model as published in the paper "Modeling geographic vaccination strategies for COVID-19 in Norway"

This repository include the code for the metapopulation model used, for the IBM code see [here](https://github.com/folkehelseinstituttet/COVID19_vaccination-IBM).

The model is calibrated using the fit_baseline_lhs.R file, then additional calibration of geographic factors to account for changes in how the metapopulations are defined when we change regional prioiritisation can be found in fit_geo_betas.R.

Once the model is calibrated, all the different vaccination scenarios can be run by running the run_analysis.R file. Code for plotting and processing of the results can be found with the IBM code.

**Note:** The data on the number of people vaccinated by day, municipality, age group and risk factor the model needs can not be shared in this repository due to potential indirect identification. For access to this data an application has to be sent to SYSVAK the Norwegian vaccination registry. For this repositry the needed files are "new_VaksinertPerKommune-seed_X-municip_code.csv" and "new_VaksinertPerKommune-seed_X-municip_code.csv" for X=1:100.

This repository depends on the metapopulation model published in [here](https://github.com/Gulfa/metapop) and Norway specific configuration [here](https://github.com/Gulfa/metapopnorge). These repositories and other requirements can be installed by running instal_reqs.R. 

