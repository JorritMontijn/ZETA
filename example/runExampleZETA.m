function runExampleZETA
	%runExampleZETA Run example ZETA-test
	%
	%This code loads data from an example LP cell and performs a ZETA-test,
	%makes a raster plot and calculates the instantaneous firing rate
	%
	%Version history:
	%1.0 - 15 June 2020
	%	Created by Jorrit Montijn

	%% load data for example cell
	rng(1,'twister'); % to match Python output
	sLoad = load('ExampleDataZETA.mat'); %loads matlab data file
	
	%some information about the neuron is stored in the sNeuron structure,
	%such as whether Kilosort2 thought it was an acceptable neuron
	sNeuron = sLoad.sNeuron;
	if sNeuron.KilosortGood == 0 || sNeuron.NonStationarity > 0.5
		error([mfilename ':BadUnit'],'This unit is non-stationary, noise-like, or contaminated');
	end
	
	%retrieve the spike times as a vector from the field in sNeuron
	vecSpikeTimes = sNeuron.SpikeTimes;
	
	%% load stimulation information
	sStim = sLoad.sStim;
	vecStimulusStartTimes = sStim.StimOnTime(:); %use (:) to ensure it's a column vector
	vecStimulusStopTimes = sStim.StimOffTime(:);
	
	%% put stimulus start and stop times together into a [T x 2] matrix
	matEventTimes = cat(2,vecStimulusStartTimes,vecStimulusStopTimes);
	
	%% run the ZETA-test with default parameters
	%if we simply want to know if the neuron responds, no hassle, we can
	%use this simple syntax with default parameters:
	dblZetaP_default = getZeta(vecSpikeTimes,matEventTimes);

	%% run the ZETA-test with specified parameters
	%however, we can also specify the parameters ourselves
	dblUseMaxDur = median(diff(vecStimulusStartTimes)); %median of trial-to-trial durations
	intResampNum = 50; %50 random resamplings should give us a good enough idea if this cell is responsive. If it's close to 0.05, we should increase this #.
	intPlot = 3;%what do we want to plot?(0=nothing, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot)
	intLatencyPeaks = 4; %how many latencies do we want? 1=ZETA, 2=-ZETA, 3=peak, 4=first crossing of peak half-height
	vecRestrictRange = [0 inf];%do we want to restrict the peak detection to for example the time during stimulus? Then put [0 1] here.
	boolVerbose = true;%displays the progress if it takes >5 seconds
	
	%then run ZETA with those parameters
	[dblZetaP,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolVerbose);
end
