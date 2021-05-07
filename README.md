# Florida-Goliath-Grouper-Trends
Estimating population trajectories for Goliath Grouper in the coastal waters of Florida (Florida Keys &amp; Dry Tortugas, Atlantic coast, and Gulf of Mexico) using Bayesian hierarchical state-space models and roving diver citizen science data.

# Data
Data on grouper abundance come from the Reef Environmental and Education Foundation (reef.org) Volunteer Fish Survey Project. Each survey notes the relative abundance, or absence, of species on a roving diver survey. We used data from all surveys conducted at dive sites where Goliath Grouper have been repeatedly sighted (at least 3 different years) from 1994 to 2020. This resulted in 17771 surveys of 130 sites across the state of Florida - with 12048 surveys of 58 sites in the Florida Keys and Dry Tortugas, 3542 surveys of 45 sites on the east coast (Biscayne NP northwards), and 389 surveys of 15 sites in the Gulf of Mexico. 

# Model
We built a Bayesian hierarchical state-space model, that estimates the underlying population trajectory over time based on the frequency of encounters and their relative abundance score. This model accounts for site-level differences in abundance (eg. sites that consistently host groupers and those that infrequently do), diver-level ability and effort, spatiotemporal clustering (eg. repeated surveys at the same sites on the same days), and detection covariates (eg. visibility, depth, etc.). The model produces estimates of annual estimated counts per survey and the underlying population index through time.

#Figures
See the figures folder for the model ouputs indicating the population index of Goliath Grouper from 1994 to 2020 for all of Florida, the Florida Keys and Dry Tortugas, and the east coast of Florida. We did not separately model the Gulf of Mexico sites due to low sample sizes.
