function [dblZetaP,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolDirectQuantile,dblJitterSize)
	%getZeta Calculates neuronal responsiveness index zeta
	%syntax: [dblZetaP,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolDirectQuantile,dblJitterSize)
	%	input:
	%	- vecSpikeTimes [S x 1]: spike times (in seconds)
	%	- vecEventTimes [T x 1]: event on times (s), or [T x 2] including event off times to calculate mean-rate difference
	%	- dblUseMaxDur: float (s), window length for calculating ZETA: ignore all spikes beyond this duration after event onset
	%								[default: minimum of all event onsets to next event onset]
	%	- intResampNum: integer, number of resamplings (default: 100)
	%	- intPlot: integer, plotting switch (0=none, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot) (default: 0)
	%	- intLatencyPeaks: integer, maximum number of latency peaks to return (1-4) (default: 2)
	%	- vecRestrictRange: temporal range within which to restrict onset/peak latencies (default: [-inf inf])
	%	- boolDirectQuantile; boolean, switch to use the empirical null-distribution rather than the
	%								Gumbel approximation (default: false) [Note: requires many resamplings!]
	%	- dblJitterSize; scalar, sets the temporal jitter window relative to dblUseMaxDur (default: 2)
	%
	%	output:
	%	- dblZetaP; p-value based on Zenith of Event-based Time-locked Anomalies
	%	- vecLatencies; different latency estimates, number determined by intLatencyPeaks. If no peaks are detected, it returns NaNs:
	%		1) Latency of ZETA
	%		2) Latency of largest z-score with inverse sign to ZETA
	%		3) Peak time of instantaneous firing rate
	%		4) Onset time of above peak, defined as the first crossing of peak half-height
	%	- sZETA; structure with fields:
	%		- dblZETA; responsiveness z-score (i.e., >2 is significant)
	%		- dblD; temporal deviation value underlying ZETA
	%		- dblP; p-value corresponding to ZETA
	%		- dblPeakT; time corresponding to ZETA
	%		- intPeakIdx; entry corresponding to ZETA
	%		- dblMeanD; Cohen's D based on mean-rate stim/base difference
	%		- dblMeanP; p-value based on mean-rate stim/base difference
	%		- vecSpikeT: timestamps of spike times (corresponding to vecD)
	%		- vecD; temporal deviation vector of data
	%		- matRandD; baseline temporal deviation matrix of jittered data
	%		- dblD_InvSign; largest peak of inverse sign to ZETA (i.e., -ZETA)
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
	%v2.10 - 3 Nov 2021
	
	%Version history:
	%0.9 - 27 June 2019
	%	Created by Jorrit Montijn
	%1.0 - 24 September 2019
	%	New procedure to determine statistical significance [by JM]
	%2.0 - 27 January 2020
	%	New peak detection procedure using multi-scale derivatives [by JM]
	%2.1 - 5 February 2020
	%	Minor changes and bug fixes [by JM]
	%2.2 - 11 February 2020
	%	Peak width, analytical ZETA correction [by JM]
	%2.3 - 26 February 2020
	%	MSD-based instantaneous spiking rates, onset latency [by JM]
	%2.4 - 11 March 2020
	%	Statistical significance using Gumbel approximation [by	JM]
	%2.5 - 27 May 2020
	%	Standardized syntax and variable names [by JM]
	%2.6 - 27 Nov 2020
	%	Improved computation time; now uses parallel bootstrapping [by JM]
	%	In case of only requesting dblZetaP, computation is now up to 10x faster
	%2.7 - 21 Jan 2021
	%	Improved computation time, using calcZeta/getSpikeT subfunctions [by JM]
	%2.8 - 23 Sept 2021
	%	Added switch to use empirical null distribution for significance calculation [by JM]
	%2.9 - 29 Oct 2021
	%	Added option to change the jitter size [by JM]
	%2.10 - 3 Nov 2021
	%	Fixed figure maximization in new matlab versions that deprecated javaframe functionality [by JM]
	
	%% prep data
	%ensure orientation
	vecSpikeTimes = vecSpikeTimes(:);
	assert(isnumeric(vecSpikeTimes),[mfilename ':WrongInputType'], 'Supplied spike time variable is not a numeric vector');
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblMeanD = nan;
	if size(matEventTimes,2) > 2
		matEventTimes = matEventTimes';
	end
	if size(matEventTimes,2) == 2
		boolStopSupplied = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = min(diff(matEventTimes(:,1)));
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 100;
	end
	
	%get intPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get intLatencyPeaks
	if ~exist('intLatencyPeaks','var') || isempty(intLatencyPeaks)
		if nargout > 1
			intLatencyPeaks = 2;
		else
			intLatencyPeaks = 0;
		end
	end
	
	%get vecRestrictRange
	if ~exist('vecRestrictRange','var') || isempty(vecRestrictRange)
		vecRestrictRange = [-inf inf];
	end

	%get boolDirectQuantile
	if ~exist('boolDirectQuantile','var') || isempty(boolDirectQuantile)
		boolDirectQuantile = false;
	end
	
	%get dblJitterSize
	if ~exist('dblJitterSize','var') || isempty(dblJitterSize)
		dblJitterSize = 2;
	end
	
	%% get zeta
	vecEventStarts = matEventTimes(:,1);
	if numel(vecEventStarts) > 1 && numel(vecSpikeTimes) > 1 && ~isempty(dblUseMaxDur) && dblUseMaxDur>0
		[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
			calcZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize);
	else
		intZETALoc = nan;
	end
	
	%% build placeholder outputs
	sRate = [];
	sZETA = [];
	vecLatencies = [];
	if isnan(intZETALoc)
		dblZetaP = 1;
		dblZETA = 0;
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate zeta');
		
		%build placeholder outputs
		if numel(vecLatencies) < intLatencyPeaks
			vecLatencies(end+1:intLatencyPeaks) = nan;
		end
		sZETA = struct;
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = 0;
		sZETA.dblP = 1;
		sZETA.dblPeakT = nan;
		sZETA.intPeakIdx = [];
		if boolStopSupplied
			sZETA.dblMeanD = 0;
			sZETA.dblMeanP = 1;
		end
		sZETA.vecSpikeT = [];
		sZETA.vecD = [];
		sZETA.matRandD = [];
		
		sZETA.dblD_InvSign = 0;
		sZETA.dblPeakT_InvSign = nan;
		sZETA.intPeakIdx_InvSign = [];
		sZETA.dblUseMaxDur = nan;
		sZETA.vecLatencyVals = [];
		return
	end
	
	%% extract real outputs
	%get location
	dblMaxDTime = vecSpikeT(intZETALoc);
	dblD = vecRealDiff(intZETALoc);
	
	%find peak of inverse sign
	[dummy,intPeakLocInvSign] = max(-sign(dblD)*vecRealDiff);
	dblMaxDTimeInvSign = vecSpikeT(intPeakLocInvSign);
	dblD_InvSign = vecRealDiff(intPeakLocInvSign);
	
	%% calculate mean-rate difference with t-test
	if boolStopSupplied && (nargout > 2 || intPlot > 1)
		vecRespBinsDur = sort(flat([matEventTimes(:,1) matEventTimes(:,2)]));
		vecR = histcounts(vecSpikeTimes,vecRespBinsDur);
		vecD = diff(vecRespBinsDur)';
		vecMu_Dur = vecR(1:2:end)./vecD(1:2:end);
		dblStart1 = min(vecRespBinsDur);
		dblFirstPreDur = dblStart1 - max([0 dblStart1 - median(vecD(2:2:end))]);
		dblR1 = sum(vecSpikeTimes > (dblStart1 - dblFirstPreDur) & vecSpikeTimes < dblStart1);
		vecMu_Pre = [dblR1 vecR(2:2:end)]./[dblFirstPreDur vecD(2:2:end)];
		
		%get metrics
		dblMeanD = mean(vecMu_Dur - vecMu_Pre) / ( (std(vecMu_Dur) + std(vecMu_Pre))/2);
		[h,dblMeanP,ci,stats]=ttest(vecMu_Dur,vecMu_Pre);
		dblMeanZ = -norminv(dblMeanP/2);
	end
	
	%% plot
	if intPlot > 1
		%plot maximally 50 traces
		intPlotIters = min([size(matRandDiff,2) 50]);
		
		%maximize figure
		figure;
		drawnow;
		try
			try
				%try new method
				h = handle(gcf);
				h.WindowState = 'maximized';
			catch
				%try old method with javaframe (deprecated as of R2021)
				sWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				drawnow;
				warning(sWarn);
			end
		catch
		end
		if intPlot > 2
			subplot(2,3,1)
			plotRaster(vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur,10000);
			xlabel('Time after event (s)');
			ylabel('Trial #');
			title('Spike raster plot');
			fixfig;
			grid off;
		end
		
		%plot
		subplot(2,3,2)
		sOpt = struct;
		sOpt.handleFig =-1;
		if dblUseMaxDur < 0.5
			dblBinSize = dblUseMaxDur/40;
		else
			dblBinSize = 0.025;
		end
		vecBins = 0:dblBinSize:dblUseMaxDur;
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,vecBins,vecEventStarts(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecMean,vecSEM);
		ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Mean spiking over trials'));
		xlabel('Time after event (s)');
		ylabel('Mean spiking rate (Hz)');
		fixfig
		
		subplot(2,3,3)
		plot(vecSpikeT,vecRealFrac)
		hold on
		plot(vecSpikeT,vecRealFracLinear,'color',[0.5 0.5 0.5]);
		title(sprintf('Real data'));
		xlabel('Time after event (s)');
		ylabel('Fractional position of spike in trial');
		fixfig
		
		subplot(2,3,4)
		cla;
		hold all
		for intOffset=1:intPlotIters
			plot(vecSpikeT,matRandDiff(:,intOffset),'Color',[0.5 0.5 0.5]);
		end
		plot(vecSpikeT,vecRealDiff,'Color',lines(1));
		scatter(dblMaxDTime,vecRealDiff(intZETALoc),'bx');
		scatter(dblMaxDTimeInvSign,vecRealDiff(intPeakLocInvSign),'b*');
		hold off
		xlabel('Time after event (s)');
		ylabel('Offset of data from linear (s)');
		if boolStopSupplied
			title(sprintf('ZETA=%.3f (p=%.3f), d(Hz)=%.3f (p=%.3f)',dblZETA,dblZetaP,dblMeanD,dblMeanP));
		else
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblZetaP));
		end
		fixfig
	end
	
	%% calculate instantaneous firing rates
	if intLatencyPeaks > 0 || nargout > 3 || intPlot > 0
		%get average of multi-scale derivatives, and rescaled to instantaneous spiking rate
		dblMeanRate = (numel(vecSpikeT)/(dblUseMaxDur*numel(vecEventStarts)));
		[vecRate,sRate] = getMultiScaleDeriv(vecSpikeT,vecRealDiff,[],[],[],intPlot,dblMeanRate,dblUseMaxDur);
	end
	
	%% calculate IFR statistics
	if ~isempty(sRate) && intLatencyPeaks > 0
		%get IFR peak
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
				vecLatencies = [dblMaxDTime dblMaxDTimeInvSign dblPeakTime dblOnset];
				vecLatencyVals = [vecRate(intZETALoc) vecRate(intPeakLocInvSign) vecRate(intPeakLoc) dblOnsetVal];
			else
				sRate.dblOnset = [nan];
				vecLatencies = [dblMaxDTime dblMaxDTimeInvSign dblPeakTime];
				vecLatencyVals = [vecRate(intZETALoc) vecRate(intPeakLocInvSign) vecRate(intPeakLoc)];
			end
			intLatencyPeaks = min([intLatencyPeaks numel(vecLatencies)]);
			vecLatencies = vecLatencies(1:intLatencyPeaks);
			vecLatencyVals = vecLatencyVals(1:intLatencyPeaks);
			if intPlot > 0
				hold on
				scatter(dblPeakTime,vecRate(intPeakLoc),'gx');
				scatter(dblMaxDTime,vecRate(intZETALoc),'bx');
				scatter(dblMaxDTimeInvSign,vecRate(intPeakLocInvSign),'b*');
				if intLatencyPeaks > 3
					scatter(dblOnset,dblOnsetVal,'rx');
					title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms,On=%.2fms',dblMaxDTime*1000,dblMaxDTimeInvSign*1000,dblPeakTime*1000,dblOnset*1000));
				else
					title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms',dblMaxDTime*1000,dblMaxDTimeInvSign*1000,dblPeakTime*1000));
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
					plot(dblMaxDTime*[1 1],vecY,'b--');
					plot(dblMaxDTimeInvSign*[1 1],vecY,'b-.');
					hold off
				end
			end
		else
			%placeholder peak data
			sRate.dblOnset = [nan];
			vecLatencies = [nan nan nan nan];
			vecLatencyVals = [nan nan nan nan];
		end
	else
		vecLatencies = [];
		vecLatencyVals = [];
	end
	
	%check number of latencies
	if numel(vecLatencies) < intLatencyPeaks
		vecLatencies(end+1:intLatencyPeaks) = nan;
		vecLatencyVals(end+1:intLatencyPeaks) = nan;
	end
	
	%% build optional output structure
	if nargout > 2
		sZETA = struct;
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = dblD;
		sZETA.dblP = dblZetaP;
		sZETA.dblPeakT = dblMaxDTime;
		sZETA.intPeakIdx = intZETALoc;
		if boolStopSupplied
			sZETA.dblMeanZ = dblMeanZ;
			sZETA.dblMeanP = dblMeanP;
		end
		sZETA.vecSpikeT = vecSpikeT;
		sZETA.vecD = vecRealDiff;
		sZETA.matRandD = matRandDiff;
		
		sZETA.dblD_InvSign = dblD_InvSign;
		sZETA.dblPeakT_InvSign = dblMaxDTimeInvSign;
		sZETA.intPeakIdx_InvSign = intPeakLocInvSign;
		sZETA.dblUseMaxDur = dblUseMaxDur;
		sZETA.vecLatencyVals = vecLatencyVals;
	end
end

