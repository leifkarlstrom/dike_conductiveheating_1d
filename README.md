# dike_conductiveheating_1d
**Conductive heating/melting model used in Biasi and Karlstrom (2021)**

Matlab codes associated with conductive heating model for Magnetogeothermometry study of Columbia River Basalt dikes. Predicts the distance over which temperatures around a dike exceed the Curie point temperature, as a function of intrusion active lifetime. 

**Instructions for use:** To conduct parameter sweeps of Currie temperature resetting as in figures S29 and S30, please run wrapper_dikeCurriePt_plots.m

diffusion_1d_dikewrapper.m could be used for single runs of the diffusion model

diffusion_1d_dike.m is the primary solver

meltfraction.m contains material models parameterized as melt fraction/temperature curves

Copyright 2021 Leif Karlstrom

