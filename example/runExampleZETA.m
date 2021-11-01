%runExampleZETA Run example ZETA-test
%
%This code loads data from an example LP cell and performs a ZETA-test,
%makes a raster plot and calculates the instantaneous firing rate
%
%Version history:
%1.0 - 15 June 2020
%	Created by Jorrit Montijn
%1.1 - 23 Sept 2021
%	Added switch to use empirical null distribution for significance calculation [by JM]
	
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
matEventTimes = cat(2,vecStimulusStartTimes,vecStimulusStopTimes); % put stimulus start and stop times together into a [T x 2] matrix

%% calculate instantaneous firing rate without performing the ZETA-test
%if we simply want to plot the neuron's response, we can use:
[vecTime,vecRate,sIFR] = getIFR(vecSpikeTimes,vecStimulusStartTimes);

%% run the ZETA-test with default parameters
%if we simply want to know if the neuron responds, no hassle, we can
%use this simple syntax with default parameters:
dblZetaP_default = getZeta(vecSpikeTimes,vecStimulusStartTimes);

%% run the ZETA-test with specified parameters
%however, we can also specify the parameters ourselves
dblUseMaxDur = median(diff(vecStimulusStartTimes)); %median of trial-to-trial durations
intResampNum = 50; %50 random resamplings should give us a good enough idea if this cell is responsive. If it's close to 0.05, we should increase this #.
intPlot = 3;%what do we want to plot?(0=nothing, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot)
intLatencyPeaks = 4; %how many latencies do we want? 1=ZETA, 2=-ZETA, 3=peak, 4=first crossing of peak half-height
vecRestrictRange = [0 inf];%do we want to restrict the peak detection to for example the time during stimulus? Then put [0 1] here.
boolDirectQuantile = false;%if true; uses the empirical null distribution rather than the Gumbel approximation. Note that in this case the accuracy of your p-value is limited by the # of resamplings

%then run ZETA with those parameters
[dblZetaP,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolDirectQuantile);

%% by popular demand: using a baseline that precedes the onset
%putting the baseline before the stimulus can be done by simply subtracting
%the baseline duration from the event times:
dblBaselineDuration = 0.5;
matEventTimesWithPrecedingBaseline = matEventTimes - dblBaselineDuration;

%then run ZETA with the new times
[dblZetaP_pb,vecLatencies_pb,sZETA_pb,sRate_pb] = getZeta(vecSpikeTimes,matEventTimesWithPrecedingBaseline,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolDirectQuantile);

%% however, the zeta function of course won't be able to tell the difference, so all timings are off by 500 ms.
%here we change the figure labels/titles (you can ignore this bit if you're not using the figure)
dblBaselineDurationMs = dblBaselineDuration*1000;
drawnow;hFig = gcf;
for intPlot=1:numel(hFig.Children)
	%adjust x-ticks
	if contains(hFig.Children(intPlot).XLabel.String,'Time ')
		set(hFig.Children(intPlot),'xticklabel',cellfun(@(x) num2str(str2double(x)-dblBaselineDuration),get(hFig.Children(intPlot),'xticklabel'),'UniformOutput',false));
	end
	%adjust timings in title
	strTitle = hFig.Children(intPlot).Title.String;
	[vecStart,vecStop]=regexp(strTitle,'[=].*?[m][s]');
	for intEntry=1:numel(vecStart)
		strOldNumber=hFig.Children(intPlot).Title.String((vecStart(intEntry)+1):(vecStop(intEntry)-2));
		strTitle = strrep(strTitle,strcat('=',strOldNumber,'ms'),strcat('=',num2str(str2double(strOldNumber)-dblBaselineDurationMs),'ms'));
	end
	hFig.Children(intPlot).Title.String = strTitle;
end

%here we adjust the times in the variables that getZeta returns
vecLatencies_pb = vecLatencies_pb - dblBaselineDuration;
sZETA_pb.vecSpikeT = sZETA_pb.vecSpikeT - dblBaselineDuration;
sRate_pb.vecT = sRate_pb.vecT - dblBaselineDuration;
sRate_pb.dblPeakTime = sRate_pb.dblPeakTime - dblBaselineDuration;
sRate_pb.dblOnset = sRate_pb.dblOnset - dblBaselineDuration;