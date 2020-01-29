# ZETA
Repository containing ZETA functions and dependencies
 
This repository contains four main functions
1) getZeta.m: Calculates the Zenith of Event-based Time-locked Anomalies (ZETA) for spike times of a single neuron.
2) getTraceZeta.m: Calculates ZETA for trace-based data, such as a calcium imaging's dF/F0 or a patch-recording's membrane voltage.
3) getMultiScaleDeriv.m: Calculates multi-scale derivatives for trace-based data, such as spike-time/z-score combinations that underlie ZETA.
4) getMultiScaleSpikeDeriv.m: Wrapper function for getMultiScaleDeriv.m when the data is spike times and event times.

Rationale for ZETA

The validity of neurophysiological studies depends on a reliable quantification of wheth-er and when a neuron responds to stimulation, be it sensory, optogenetically or other-wise. However, current statistical analysis methods to determine a neuron’s respon-siveness are often labour-intensive and/or only detect classically mean-rate modulated cells. This problem is becoming more acute with the recent advent of techniques that yield increasingly large numbers of cells, such as multi-plane calcium imaging and Neu-ropixels recordings. Moreover, using peristimulus time histograms (PSTHs) to determine a neuron’s responsiveness still requires an a priori selection of an appropriate bin size, which can be heterogeneous over neurons. The procedure presented here, ZETA (Zen-ith of Event-based Time-locked Anomalies), consistently and robustly outperforms common approaches for quantifying neuronal responsiveness, in the sense that it includes more cells at a similar false-positive rate. This effect holds in both artificially generated benchmarking datasets with known ground truths, as well as in experimentally recorded electrophysiological data (see /ZETA/benchmarks/). Our procedure automatically includes otherwise undetectable non-trivially modulated neurons from a variety of brain regions and recording techniques, including Neuropixels recordings in mouse visual cortex and subcortical area SC (superior colliculus), as well as retinal ganglion cell spiking responses to light flashes recorded with multielectrode arrays, and GCaMP6f imaging in primary visual cortex with natural movies and drifting gratings. 

Finally, ZETA’s timescale-, parameter- and binning-free nature allowed us to implement a ZETA-derived algorithm (using multi-scale derivatives) to calculate peak onset and offset latencies in neuronal spike trains with theoretically unlimited temporal resolution. 

Please send any questions or comments to j.montijn at nin.knaw.nl.
