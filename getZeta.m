function [dblZETA,vecLatencies,sZETA,sMSD] = getZeta(vecSpikeTimes,varEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,boolVerbose)
	%getZeta Calculates neuronal responsiveness index zeta
	%syntax: [dblZETA,vecLatencies,sZETA,sMSD] = getZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,intPlot,boolVerbose)
	%	input:
	%	- vecSpikeTimes [S x 1]: spike times (in seconds)
	%	- vecEventTimes [T x 1]: event on times (s), or [T x 2] including event off times to calculate mean-rate difference
	%	- dblUseMaxDur: float (s), window length for calculating ZETA: ignore all spikes beyond this duration after event onset
	%								[default: median of event onset to event onset]
	%	- intResampNum: integer, number of resamplings (default: 50)
	%	- intPlot: integer, plotting switch (0=none, 1=traces only, 2=raster plot as well) (default: 0)
	%	- intLatencyPeaks: integer, maximum number of latency peaks to return (default: 3)
	%	- boolVerbose: boolean, switch to print messages
	%
	%	output:
	%	- dblZETA; Zenith of Event-based Time-locked Anomalies: FDR-corrected responsiveness z-score (i.e., >2 is significant)
	%	- vecLatencies; different latency estimates, number determined by intLatencyPeaks. If no peaks are detected, it returns NaNs:
	%		1) Latency of ZETA
	%		2) Latency of largest z-score with inverse sign to ZETA
	%		3) Peak time of multi-scale derivatives
	%		4...N) Onsets of largest sustained peaks
	%	- sZETA; structure with fields:
	%		- dblZ; uncorrected peak z-score 
	%		- dblP; p-value corresponding to zeta
	%		- dblPeakT; time corresponding to ZETA
	%		- dblHzD; Cohen's D based on mean-rate stim/base difference
	%		- dblHzP; p-value based on mean-rate stim/base difference
	%		- vecSpikeT: timestamps of spike times (corresponding to vecZ)
	%		- vecZ; z-score for all time points corresponding to vecSpikeT
	%		- vecRealDiff: real offset of spikes relative to uniform rate
	%		- matRandDiff; matrix of shuffled runs with offset to uniform
	%	- sMSD; structure with fields: (only if intLatencyPeaks > 0)
	%		- vecMSD; Multi-scale derivative provides high-resolution time-locked activation (like a PSTH)
	%		- vecScale; timescales used to calculate derivatives
	%		- matSmoothMSD; smoothed multi-scale derivatives matrix
	%		- matMSD; raw multi-scale derivatives matrix
	%		Data on the peak:
	%		- dblPeakTime; time of peak (in seconds)
	%		- dblPeakWidth; duration of peak (in seconds)
	%		- vecPeakStartStop; start and stop time of peak (in seconds)
	%		- intPeakLoc; spike index of peak (corresponding to sZETA.vecSpikeT)
	%		- vecPeakStartStopIdx; spike indices of peak start/stop (corresponding to sZETA.vecSpikeT)
	%		Additionally, it will return peak latencies using findsustainedpeaks.m:
	%		- vecTime: latencies for activation
	%		- vecIdx:  indices corresponding to peak locations in vecSpikeT/vecZ entries
	%		- vecDuration: peak durations
	%		- vecEnergy: peak energies
	%
	%Version history:
	%0.9 - June 27 2019
	%	Created by Jorrit Montijn
	%1.0 - September 24 2019
	%	New procedure to determine statistical significance [by JM]
	%2.0 - January 27 2020
	%	New peak detection procedure using multi-scale derivatives [by JM]
	%2.1 - February 5 2020
	%	Minor changes and bug fixes [by JM]
	%2.2 - February 11 2020
	%	Peak width, analytical ZETA correction [by JM]
	
	%% prep data
	%ensure orientation
	vecSpikeTimes = vecSpikeTimes(:);
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblHzD = nan;
	if size(varEventTimes,2) > 2
		varEventTimes = varEventTimes';
	end
	if size(varEventTimes,2) == 2
		boolStopSupplied = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(varEventTimes(:,1)));
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 50;
	end
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get boolPlot
	if ~exist('intLatencyPeaks','var') || isempty(intLatencyPeaks)
		intLatencyPeaks = 3;
	end
	
	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = false;
	end
	
	%% build onset/offset vectors
	vecEventStarts = varEventTimes(:,1);
	vecEventStops = varEventTimes(:,2);
	
	%% prepare interpolation points
	%pre-allocate
	intMaxRep = size(varEventTimes,1);
	cellSpikeTimesPerTrial = cell(intMaxRep,1);
	
	%go through trials to build spike time vector
	for intEvent=1:intMaxRep
		%get times
		dblStartT = vecEventStarts(intEvent,1);
		dblStopT = dblStartT + dblUseMaxDur;
		
		% build trial assignment
		cellSpikeTimesPerTrial{intEvent} = vecSpikeTimes(vecSpikeTimes < dblStopT & vecSpikeTimes > dblStartT) - dblStartT;
	end
	
	%get spikes in fold
	vecSpikeT = sort(cell2vec(cellSpikeTimesPerTrial),'ascend');
	intSpikes = numel(vecSpikeT);
	
	%% run normal
	%get data
	[vecRealDiff,vecRealFrac,vecRealFracLinear] = ...
		getTempOffset(vecSpikeT,vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);

	%% run bootstraps2
	hTic = tic;
	matRandDiff = nan(intSpikes,intResampNum);
	for intResampling=1:intResampNum
		%% msg
		if boolVerbose && toc(hTic) > 5
			fprintf('Now at resampling %d/%d\n',intResampling,intResampNum);
			hTic = tic;
		end
		%% get random subsample
		vecStimUseOnTime = vecEventStarts(:,1) + 2*dblUseMaxDur*(rand(size(vecEventStarts(:,1)))-0.5);
		
		%get temp offset
		[vecRandDiff,vecRandFrac,vecRandFracLinear] = ...
			getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur);
	
		%assign data
		matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
	end

	%% calculate measure of effect size (for equal n, d' equals Cohen's d)
	%define plots
	dblRandSd = nanstd(matRandDiff(:));
	vecZ = vecRealDiff./dblRandSd;
	if numel(vecZ) < 3
		dblZETA = 0;
		sZETA = struct;
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate zeta');
		return
	end
	
	%find highest peak and retrieve value
	[dummy,intPeakLoc]= max(abs(vecZ));
	dblMaxZTime = vecSpikeT(intPeakLoc);
	dblZ = vecZ(intPeakLoc);
	dblZETA = (sqrt(2)/sqrt(pi)) + dblZ*(1 - (2/pi)); %apply correction factor for half-normal
	dblP=1-(normcdf(abs(dblZETA))-normcdf(-abs(dblZETA)));
	%find peak of inverse sign
	[dummy,intPeakLocInvSign] = max(-sign(dblZ)*vecZ);
	dblMaxZTimeInvSign = vecSpikeT(intPeakLocInvSign);
	dblZ_InvSign = vecZ(intPeakLocInvSign);
	
	if boolStopSupplied
		%% calculate mean-rate difference
		%pre-allocate
		vecStimHz = zeros(intMaxRep,1);
		vecBaseHz = zeros(intMaxRep,1);
		dblMedianBaseDur = median(vecEventStarts(2:end) - vecEventStops(1:(end-1)));
		
		%go through trials to build spike time vector
		for intEvent=1:intMaxRep
			%get times
			dblStartT = vecEventStarts(intEvent,1);
			dblStopT = dblStartT + dblUseMaxDur;
			dblPreT = dblStartT - dblMedianBaseDur;
			
			% build trial assignment
			vecStimHz(intEvent) = sum(vecSpikeTimes < dblStopT & vecSpikeTimes > dblStartT)/(dblStopT - dblStartT);
			vecBaseHz(intEvent) = sum(vecSpikeTimes < dblStartT & vecSpikeTimes > dblPreT)/dblMedianBaseDur;
		end
		
		%get metrics
		dblHzD = mean(vecStimHz - vecBaseHz) / ( (std(vecStimHz) + std(vecBaseHz))/2);
		[h,dblHzP]=ttest(vecStimHz,vecBaseHz);
	end
	
	%% plot
	if intPlot
		%plot maximally 50 traces
		intPlotIters = min([size(matRandDiff,2) 50]); 
		
		%make maximized figure
		figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		if intPlot == 2
			subplot(2,3,1)
			plotRaster(vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);
			xlabel('Time from event (s)');
			ylabel('Trial #');
			title('Spike raster plot');
			fixfig;
		end
		
		%plot
		subplot(2,3,2)
		sOpt = struct;
		sOpt.handleFig =-1;
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,0:0.1:dblUseMaxDur,vecEventStarts(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecMean,vecSEM);
		ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Mean spiking over trials'));
		xlabel('Time from event (s)');
		ylabel('Mean spiking rate (Hz)');
		fixfig
		
		subplot(2,3,3)
		plot(vecSpikeT,vecRealFrac)
		hold on
		plot(vecSpikeT,vecRealFracLinear,'color',[0.5 0.5 0.5]);
		title(sprintf('Real data'));
		xlabel('Time from event (s)');
		ylabel('Fractional position of spike in trial');
		fixfig
		
		subplot(2,3,4)
		cla;
		hold all
		for intOffset=1:intPlotIters
			plot(vecSpikeT,matRandDiff(:,intOffset),'Color',[0.5 0.5 0.5]);
		end
		plot(vecSpikeT,vecRealDiff,'Color',lines(1));
		scatter(dblMaxZTime,vecRealDiff(intPeakLoc),'bx');
		scatter(dblMaxZTimeInvSign,vecRealDiff(intPeakLocInvSign),'b*');
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		if boolStopSupplied
			title(sprintf('ZETA=%.3f (p=%.3f), d(Hz)=%.3f (p=%.3f)',dblZETA,dblP,dblHzD,dblHzP));
		else
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblP));
		end
		fixfig
		
	end
	
	%% calculate MSD if significant
	if intLatencyPeaks > 0
		[vecMSD,sMSD] = getMultiScaleDeriv(vecSpikeT,vecRealDiff,[],[],[],intPlot);
	else
		sMSD = [];
	end
	
	%% calculate MSD statistics
	if ~isempty(sMSD) && intLatencyPeaks > 0
		%get sustained peak onset
		[vecPeakTime,vecPeakIdx,vecPeakDuration,vecPeakEnergy,vecSlope] = findsustainedpeaks(vecMSD,vecSpikeT,intLatencyPeaks);
		%get MSD most prominent peak time
		[vecValsPos,vecLocsPos,vecsWidthPos,vecPromsPos]=findpeaks(vecMSD);
		[dblMaxPosVal,intPosIdxMSD] = max(vecValsPos);
		[vecValsNeg,vecLocsNeg,vecsWidthNeg,vecPromsNeg]=findpeaks(-vecMSD);
		[dblMaxNegVal,intNegIdxMSD] = max(vecValsNeg);
		if abs(dblMaxPosVal) > abs(dblMaxNegVal)
			intIdxMSD = intPosIdxMSD;
			intPeakLocMSD = vecLocsPos(intIdxMSD);
			dblPeakPromMSD = vecPromsPos(intIdxMSD);
			dblCutOff = vecMSD(intPeakLocMSD) - dblPeakPromMSD/2;
			indPeakMembers = vecMSD > dblCutOff;
		else
			intIdxMSD = intNegIdxMSD;
			intPeakLocMSD = vecLocsNeg(intIdxMSD);
			dblPeakPromMSD = vecPromsNeg(intIdxMSD);
			
			dblCutOff = vecMSD(intPeakLocMSD) - dblPeakPromMSD/2;
			indPeakMembers = vecMSD < dblCutOff;
		end
		%get potential starts/stops
		vecPeakStarts = find(diff(indPeakMembers)==1);
		vecPeakStops = find(diff(indPeakMembers)==-1);
		if indPeakMembers(1) == 1,vecPeakStarts = [1 vecPeakStarts(:)'];end
		if indPeakMembers(end) == 1,vecPeakStops = [vecPeakStops(:)' numel(indPeakMembers)];end
		%find closest points
		intPeakStart = intPeakLocMSD-min(intPeakLocMSD - vecPeakStarts(vecPeakStarts<intPeakLocMSD));
		intPeakStop = intPeakLocMSD+min(vecPeakStops(vecPeakStops>intPeakLocMSD) - intPeakLocMSD);
		dblPeakStartT = vecSpikeT(intPeakStart);
		dblPeakStopT = vecSpikeT(intPeakStop);
		dblPeakWidthMSD = dblPeakStopT-dblPeakStartT;
		dblPeakTimeMSD = vecSpikeT(intPeakLocMSD);
		%assign peak data
		sMSD.dblPeakTime = dblPeakTimeMSD;
		sMSD.dblPeakWidth = dblPeakWidthMSD;
		sMSD.vecPeakStartStop = [dblPeakStartT dblPeakStopT];
		sMSD.intPeakLoc = intPeakLocMSD;
		sMSD.vecPeakStartStopIdx = [intPeakStart intPeakStop];
		%assign array data
		sMSD.vecTime = vecPeakTime;
		sMSD.vecIdx = vecPeakIdx;
		sMSD.vecDuration = vecPeakDuration;
		sMSD.vecEnergy = vecPeakEnergy;
		vecLatencies = [dblMaxZTime dblMaxZTimeInvSign dblPeakTimeMSD vecPeakTime(:)'];
		vecLatencies = vecLatencies(1:intLatencyPeaks);
		if intPlot > 0
			hold on
			scatter(vecPeakTime(1),vecMSD(vecPeakIdx(1)),'rx');
			scatter(dblPeakTimeMSD,vecMSD(intPeakLocMSD),'gx');
			scatter(dblMaxZTime,vecMSD(intPeakLoc),'bx');
			scatter(dblMaxZTimeInvSign,vecMSD(intPeakLocInvSign),'b*');
			hold off
			title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms,Sh=%.0fms',dblMaxZTime*1000,dblMaxZTimeInvSign*1000,dblPeakTimeMSD*1000,vecPeakTime(1)*1000));
			fixfig;
		

			vecHandles = get(gcf,'children');
			ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
			axes(ptrFirstSubplot);
			vecY = get(gca,'ylim');
			hold on;
			plot(vecPeakTime(1)*[1 1],vecY,'r--');
			plot(dblPeakTimeMSD*[1 1],vecY,'g--');
			plot(dblMaxZTime*[1 1],vecY,'b--');
			plot(dblMaxZTimeInvSign*[1 1],vecY,'b-.');
			
			hold off
		end
	else
		vecLatencies = [];
	end
			
	%check number of latencies
	if numel(vecLatencies) < intLatencyPeaks
		vecLatencies(end+1:intLatencyPeaks) = nan;
	end

	%% build optional output structure
	if nargin > 1
		sZETA = struct;
		sZETA.dblZ = dblZ;
		sZETA.dblP = dblP;
		sZETA.dblPeakT = dblMaxZTime;
		sZETA.intPeakIdx = intPeakLoc;
		if boolStopSupplied
			sZETA.dblHzD = dblHzD;
			sZETA.dblHzP = dblHzP;
		end
		sZETA.vecSpikeT = vecSpikeT;
		sZETA.vecZ = vecZ;
		sZETA.vecRealDiff = vecRealDiff;
		sZETA.matRandDiff = matRandDiff;
		
		sZETA.dblZ_InvSign = dblZ_InvSign;
		sZETA.intPeakIdx_InvSign = intPeakLocInvSign;
		sZETA.dblPeakT_InvSign = dblMaxZTimeInvSign;
	end
end

