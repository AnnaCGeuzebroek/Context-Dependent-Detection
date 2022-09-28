# Context dependent target detection within a continuous context
> Geuzebroek AC, Craddock H, O’Connell RG, &amp; Kelly SP (2022). Balancing true and false detection of intermittent sensory targets by adjusting the inputs to the evidence accumulation process. doi: https://biorxiv.org/cgi/content/short/2022.09.01.505650v1

Code as described in Geuzebroek et al. (2022), allowing the extraction of behavioural and electrophysiological data (contextAnalysis) as well as behavioural and neurally-informed modelling of the data (contextModelling). All preprocessed data is available at DOI:[/10.17605/OSF.IO/YJVKU](https://osf.io/yjvku/?view_only=7ed5aee5d09a4d5ca13de1ba169b0588), but raw data will be made available at request. 

For an excellent review on neurally-informed modelling see:

- O’Connell, R.G., Shadlen, M.N., Wong-Lin, K., & Kelly, S.P. (2018). Bridging Neural and Computational Viewpoints on Perceptual Decision-Making. *Trends in neurosciences*. [link](https://www.sciencedirect.com/science/article/pii/S0166223618301668)

## Presentation code
Code is provided to run a continuous detection version of the random dot motion task were participants continuously monitor a cloud of white, randomly moving dots for intermittent targets defined as a step change from random to coherent upwards dot motion. During periods of incoherent motion (0% motion coherence), all dots were randomly displaced to a new location throughout the patch on each frame. Coherent dot motion was accomplished by displaying a certain percentage of randomly selected dots in a direction relative to their previous location within each frame. Code is depending utilizing Psychtoolbox functions. 

## contextAnalysis
Preprocessed data can be found [here](https://osf.io/yjvku/?view_only=7ed5aee5d09a4d5ca13de1ba169b0588). How this data was preprocessed, how behavioura data and EEG epoch were extracted are all described in details [here](https://biorxiv.org/cgi/content/short/2022.09.01.505650v1). Raw data will be made available when requested and the preprocessed data can than be reproduced. *contextAnalysis.m* is a file that uses *dataAnalysis.m* which is available at [Neurally-Informed-Modelling](https://github.com/AnnaCGeuzebroek/Neurally-Informed-Modelling), while *dataAnalysis.m* allows you to upload your own data and apply all steps, here it skip the preprocessing steps straight to the waveform plotting and statistical analysis. 

## contextModelling
Details on how the modelling was performed is described (here)[https://biorxiv.org/cgi/content/short/2022.09.01.505650v1]. *contextModelling.m* uses *dataModelling.m* avaiable at [Neurally-Informed-Modelling](https://github.com/AnnaCGeuzebroek/Neurally-Informed-Modelling). Similar to *contextAnalysis.m* it is designed to use the tmpBehavioural data in the [preprocessed data](https://osf.io/yjvku/?view_only=7ed5aee5d09a4d5ca13de1ba169b0588). Pre-calculated model fits are also made avaible there. This modelling uses an urgency signal, determining the motor preperation through the inter-target-interval, this was for this current paper calculated in a seperated code written by Simon P. Kelly *makeUrgencyFn.m*, see Dependables and the Urgency.mat. In future iterations in Neurally-Informed-Modelling](https://github.com/AnnaCGeuzebroek/Neurally-Informed-Modelling) it will be possible to extract that in *dataAnalysis.m*. Similar, in this iteration of the code, different code is used to add the urgency signal, this might be changed is Beta as motor preperation signals has been further valided in several other tasks (see *applyNIModelling.m* and *NIDecisionModels.m*)

##
<sup>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. If you use the Software for your own research, cite the paper.</sup>

<sup>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</sup>

Anna Geuzebroek and Simon Kelly, 2022
anna.geuzebroek@ucd.ie / simon.kelly@ucd.ie
