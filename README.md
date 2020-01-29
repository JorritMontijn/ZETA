# TuningMetrics
 Repository containing several tuning curve metrics.
 
The metrics included here address the following issues: 
1) Detection of cells that show a non-trivial modulation by a stimulus, in the sense that their temporal pattern of activity may be reliable, but their average firing rate is no different during the presence or absence of a stimulus. (getZeta.m)
2) A parameter-free, information-based metric of stimulus tuning that is more sensitive than common metrics such as (1 – circular variance), or the orientation selectivity index. (getDeltaPrime.m)
3) A parameter-free smoothness-based metric for quantifying tuning bandwidth/sparseness. (getTuningRho.m)

Rationale for ZETA

The validity of neurophysiological studies depends on a reliable quantification of whether and when a neuron responds to stimulation, be it sensory, optogenetically or otherwise. However, current statistical analysis methods to determine a neuron’s responsiveness are often labour-intensive and/or only detect classically mean-rate modulated cells. This problem is becoming more acute with the recent advent of techniques that yield increasingly large numbers of cells, such as multi-plane calcium imaging and Neuropixels recordings. Here, we present a procedure, ZETA (Zenith of Event-based Time-locked Anomalies), that consistently and robustly outperforms common approaches, in the sense that it includes more cells at a similar false-positive rate. Moreover, our procedure automatically includes otherwise undetectable non-trivial responses from a variety of brain regions and recording techniques. This includes, among others, retinal ganglion cell spiking responses to light flashes and calcium imaging in primary visual cortex with natural movies. 


NOTE: Zeta is benchmarked, but these metrics are still work in progress, so all code is only provided as-is. Caveat emptor.
