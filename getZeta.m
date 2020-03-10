function [dblZETA,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,varEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolVerbose)
	%getZeta Calculates neuronal responsiveness index zeta
	%syntax: [dblZETA,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,boolVerbose)
	%	input:
	%	- vecSpikeTimes [S x 1]: spike times (in seconds)
	%	- vecEventTimes [T x 1]: event on times (s), or [T x 2] including event off times to calculate mean-rate difference
	%	- dblUseMaxDur: float (s), window length for calculating ZETA: ignore all spikes beyond this duration after event onset
	%								[default: median of event onset to event onset]
	%	- intResampNum: integer, number of resamplings (default: 50)
	%	- intPlot: integer, plotting switch (0=none, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot) (default: 0)
	%	- intLatencyPeaks: integer, maximum number of latency peaks to return (1-4) (default: 4)
	%	- vecRestrictRange: temporal range within which to restrict onset/peak latencies (default: [-inf inf])
	%	- boolVerbose: boolean, switch to print progress messages (default: false)
	%
	%	output:
	%	- dblZETA; Zenith of Event-based Time-locked Anomalies: FDR-corrected responsiveness z-score (i.e., >2 is significant)
	%	- vecLatencies; different latency estimates, number determined by intLatencyPeaks. If no peaks are detected, it returns NaNs:
	%		1) Latency of ZETA
	%		2) Latency of largest z-score with inverse sign to ZETA
	%		3) Peak time of instantaneous firing rate
	%		4) Onset time of above peak, defined as the first crossing of peak half-height
	%	- sZETA; structure with fields:
	%		- dblZ; uncorrected peak z-score
	%		- dblP; p-value corresponding to zeta
	%		- dblPeakT; time corresponding to ZETA
	%		- intPeakIdx; entry corresponding to ZETA
	%		- dblMeanD; Cohen's D based on mean-rate stim/base difference
	%		- dblMeanP; p-value based on mean-rate stim/base difference
	%		- vecSpikeT: timestamps of spike times (corresponding to vecZ)
	%		- vecZ; z-score for all time points corresponding to vecSpikeT
	%		- vecRealDiff: real offset of spikes relative to uniform rate
	%		- matRandDiff; matrix of shuffled runs with offset to uniform
	%		- dblZ_InvSign; largest peak z-score of inverse sign to ZETA
	%		- dblPeakT_InvSign; time corresponding to -ZETA
	%		- intPeakIdx_InvSign; entry corresponding to -ZETA
	%		- dblUseMaxDur; window length used to calculate ZETA
	%	- sRate; structure with fields: (only if intLatencyPeaks > 0)
	%		- vecRate; instantaneous spiking rates (like a PSTH)
	%		- vecT; time-points corresponding to vecRate (same as sZETA.vecSpikeT)
	%		- vecM; Mean of multi-scale derivatives
	%		- vecScale; timescales used to calculate derivatives
	%		- matMSD; multi-scale derivatives matrix
	%		- vecV; values on which vecRate is calculated (same as sZETA.vecZ)
	%		Data on the peak:
	%		- dblPeakTime; time of peak (in seconds)
	%		- dblPeakWidth; duration of peak (in seconds)
	%		- vecPeakStartStop; start and stop time of peak (in seconds)
	%		- intPeakLoc; spike index of peak (corresponding to sZETA.vecSpikeT)
	%		- vecPeakStartStopIdx; spike indices of peak start/stop (corresponding to sZETA.vecSpikeT)
	%		Additionally, it will return peak onset latency (first crossing of peak half-height) using getOnset.m:
	%		- dblOnset: latency for peak onset
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
	%2.3 - February 26 2020
	%	MSD-based instantaneous spiking rates, onset latency [by JM]
	
	%% prep data
	%ensure orientation
	vecSpikeTimes = vecSpikeTimes(:);
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblMeanD = nan;
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
	
	%get intLatencyPeaks
	if ~exist('intLatencyPeaks','var') || isempty(intLatencyPeaks)
		intLatencyPeaks = 4;
	end
	
	%get boolPlot
	if ~exist('vecRestrictRange','var') || isempty(vecRestrictRange)
		vecRestrictRange = [-inf inf];
	end

	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = false;
	end
	
	%% build onset/offset vectors
	vecEventStarts = varEventTimes(:,1);
	
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
	vecSpikeT = [0;sort(cell2vec(cellSpikeTimesPerTrial),'ascend');dblUseMaxDur];
	intSpikes = numel(vecSpikeT);
	
	%% run normal
	%get data
	[vecRealDiff,vecRealFrac,vecRealFracLinear] = ...
		getTempOffset(vecSpikeT,vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);
	
	%% run bootstraps
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
		sZETA = [];
		vecLatencies = [];
		sRate = [];
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate zeta');
		
		%build placeholder outputs
		if numel(vecLatencies) < intLatencyPeaks
			vecLatencies(end+1:intLatencyPeaks) = nan;
		end
		return
	end
	
	%find highest peak and retrieve value
	[dummy,intZETALoc]= max(abs(vecZ));
	dblMaxZTime = vecSpikeT(intZETALoc);
	dblZ = vecZ(intZETALoc);
	dblZETA = sign(dblZ)*((sqrt(2)/sqrt(pi)) + abs(dblZ)*(1 - (2/pi))); %apply correction factor for half-normal
	dblP=1-(normcdf(abs(dblZETA))-normcdf(-abs(dblZETA)));
	%find peak of inverse sign
	[dummy,intPeakLocInvSign] = max(-sign(dblZ)*vecZ);
	dblMaxZTimeInvSign = vecSpikeT(intPeakLocInvSign);
	dblZ_InvSign = vecZ(intPeakLocInvSign);
	
	if boolStopSupplied
		%% calculate mean-rate difference
		%pre-allocate
		vecEventStops = varEventTimes(:,2);
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
		dblMeanD = mean(vecStimHz - vecBaseHz) / ( (std(vecStimHz) + std(vecBaseHz))/2);
		[h,dblMeanP]=ttest(vecStimHz,vecBaseHz);
	end
	
	%% plot
	if intPlot > 1
		%plot maximally 50 traces
		intPlotIters = min([size(matRandDiff,2) 50]);
		
		%make maximized figure
		figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		if intPlot > 2
			subplot(2,3,1)
			plotRaster(vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur,10000);
			xlabel('Time from event (s)');
			ylabel('Trial #');
			title('Spike raster plot');
			fixfig;
			grid off;
		end
		
		%plot
		subplot(2,3,2)
		sOpt = struct;
		sOpt.handleFig =-1;
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,0:0.025:dblUseMaxDur,vecEventStarts(:,1),sOpt);
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
		scatter(dblMaxZTime,vecRealDiff(intZETALoc),'bx');
		scatter(dblMaxZTimeInvSign,vecRealDiff(intPeakLocInvSign),'b*');
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		if boolStopSupplied
			title(sprintf('ZETA=%.3f (p=%.3f), d(Hz)=%.3f (p=%.3f)',dblZETA,dblP,dblMeanD,dblMeanP));
		else
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblP));
		end
		fixfig
	end
	
	%% calculate MSD if significant
	if intLatencyPeaks > 0
		%get average of multi-scale derivatives, and rescaled to instantaneous spiking rate
		dblMeanRate = (intSpikes/(dblUseMaxDur*intMaxRep));
		[vecRate,sRate] = getMultiScaleDeriv(vecSpikeT,vecRealDiff,[],[],[],intPlot,dblMeanRate,dblUseMaxDur);
	else
		sRate = [];
	end
	
	%% calculate MSD statistics
	if ~isempty(sRate) && intLatencyPeaks > 0
		%get MSD peak
		[dblPeakRate,dblPeakTime,dblPeakWidth,vecPeakStartStop,intPeakLoc,vecPeakStartStopIdx] = getPeak(vecRate,vecSpikeT,vecRestrictRange);
		sRate.dblPeakRate = dblPeakRate;
		sRate.dblPeakTime = dblPeakTime;
		sRate.dblPeakWidth = dblPeakWidth;
		sRate.vecPeakStartStop = vecPeakStartStop;
		sRate.intPeakLoc = intPeakLoc;
		sRate.vecPeakStartStopIdx = vecPeakStartStopIdx;
		
		
		if ~isnan(dblPeakTime)
			%assign array data
			if intLatencyPeaks > 3
				%get onset
				[dblOnset,dblOnsetVal] = getOnset(vecRate,vecSpikeT,dblPeakTime,vecRestrictRange);
				sRate.dblOnset = dblOnset;
				vecLatencies = [dblMaxZTime dblMaxZTimeInvSign dblPeakTime dblOnset];
			else
				sRate.dblOnset = [];
				vecLatencies = [dblMaxZTime dblMaxZTimeInvSign dblPeakTime];
			end
			vecLatencies = vecLatencies(1:intLatencyPeaks);
			if intPlot > 0
				hold on
				scatter(dblPeakTime,vecRate(intPeakLoc),'gx');
				scatter(dblMaxZTime,vecRate(intZETALoc),'bx');
				scatter(dblMaxZTimeInvSign,vecRate(intPeakLocInvSign),'b*');
				if intLatencyPeaks > 3
					scatter(dblOnset,dblOnsetVal,'rx');
					title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms,On=%.0fms',dblMaxZTime*1000,dblMaxZTimeInvSign*1000,dblPeakTime*1000,dblOnset*1000));
				else
					title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms',dblMaxZTime*1000,dblMaxZTimeInvSign*1000,dblPeakTime*1000));
				end
				hold off
				fixfig;
				
				if intPlot > 3
					vecHandles = get(gcf,'children');
					ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
					axes(ptrFirstSubplot);
					vecY = get(gca,'ylim');
					hold on;
					if intLatencyPeaks > 3,plot(dblOnset*[1 1],vecY,'r--');end
					plot(dblPeakTime*[1 1],vecY,'g--');
					plot(dblMaxZTime*[1 1],vecY,'b--');
					plot(dblMaxZTimeInvSign*[1 1],vecY,'b-.');
					hold off
				end
			end
		else
			%placeholder peak data
			sRate.dblOnset = [];
			vecLatencies = [nan nan nan nan];
		end
	else
		vecLatencies = [];
	end
	
	%check number of latencies
	if numel(vecLatencies) < intLatencyPeaks
		vecLatencies(end+1:intLatencyPeaks) = nan;
	end
	
	%% build optional output structure
	if nargout > 2
		sZETA = struct;
		sZETA.dblZ = dblZ;
		sZETA.dblP = dblP;
		sZETA.dblPeakT = dblMaxZTime;
		sZETA.intPeakIdx = intZETALoc;
		if boolStopSupplied
			sZETA.dblMeanD = dblMeanD;
			sZETA.dblMeanP = dblMeanP;
		end
		sZETA.vecSpikeT = vecSpikeT;
		sZETA.vecZ = vecZ;
		sZETA.vecRealDiff = vecRealDiff;
		sZETA.matRandDiff = matRandDiff;
		
		sZETA.dblZ_InvSign = dblZ_InvSign;
		sZETA.dblPeakT_InvSign = dblMaxZTimeInvSign;
		sZETA.intPeakIdx_InvSign = intPeakLocInvSign;
		sZETA.dblUseMaxDur = dblUseMaxDur;
	end
end

