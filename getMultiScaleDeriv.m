function [vecRate,sMSD] = getMultiScaleDeriv(vecT,vecV,intSmoothSd,dblMinScale,dblBase,intPlot,dblMeanRate,dblUseMaxDur,boolUseParallel)
	%getMultiScaleDeriv Returns multi-scale derivative. Syntax:
	%   [vecRate,sMSD] = getMultiScaleDeriv(vecT,vecV,intSmoothSd,dblMinScale,dblBase,intPlot,dblMeanRate,dblUseMaxDur,boolUseParallel)
	%Required input:
	%	- vecT [N x 1]: timestamps (e.g., spike times)
	%	- vecV [N x 1]: values (e.g., z-scores)
	%
	%Optional inputs:
	%	- intSmoothSd: Gaussian SD of smoothing kernel (in # of samples) [default: 0]
	%	- dblMinScale: minimum derivative scale in seconds [default: 1/1000]
	%	- dblBase: base for exponential scale step size [default: 1.5]
	%	- intPlot: integer, plotting switch (0=none, 1=plot rates, 2=subplot 5&6 of [2 3]) [default: 0]. 
	%						If set to 1, it will plot into the current axes if empty, or create a new figure if ~isempty(get(gca,'Children'))
	%	- dblMeanRate: mean spiking rate to normalize vecRate (optional)
	%	- dblUseMaxDur: trial duration to normalize vecRate (optional)
	%	- boolUseParallel: use parallel processing (optional) [default: 0; often decreases performance, so be cautious!]
	%
	%Outputs:
	%	- vecRate; Instantaneous spiking rate
	%	- sMSD; structure with fields:
	%		- vecRate; instantaneous spiking rates (like a PSTH)
	%		- vecT; time-points corresponding to vecRate (same as input vecT)
	%		- vecM; Mean of multi-scale derivatives
	%		- vecScale; timescales used to calculate derivatives
	%		- matMSD; multi-scale derivatives matrix
	%		- vecV; values on which vecRate is calculated (same as input vecV)
	%
	%Version history:
	%1.0 - January 24 2020
	%	Created by Jorrit Montijn - split from previous getMultiScaleDeriv.
	%1.1 - February 26 2020
	%	Added instantaneous spiking rate rescaling [by JM]
	%1.1.1 - January 10 2022
	%	Changed plotting behavior to create new figure when intPlot==1 if gca is not empty [by JM]
	
	%% set default values
	if ~exist('intSmoothSd','var') || isempty(intSmoothSd)
		intSmoothSd = 0;
	end
	if ~exist('dblBase','var') || isempty(dblBase)
		dblBase = 1.5;
	end
	if ~exist('dblMinScale','var') || isempty(dblMinScale)
		dblMinScale = round(log(1/1000) / log(dblBase));
	end
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	if ~exist('dblMeanRate','var') || isempty(dblMeanRate)
		dblMeanRate = 1;
		strLabelY = 'Time-locked activation (a.u.)';
	else
		strLabelY = 'Spiking rate (Hz)';
	end
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = range(vecT);
	end
	if ~exist('boolUseParallel','var') || isempty(boolUseParallel)
		boolUseParallel = false;
	end
	
	%% reorder just in case
	[vecT,vecReorder] = sort(vecT(:),'ascend');
	vecV = vecV(vecReorder);
	vecV = vecV(:);
	
	%% prepare data
	dblMaxScale = log(max(vecT)/10) / log(dblBase);
	intN = numel(vecT);
	
	%% get multi-scale derivative
	vecExp = dblMinScale:dblMaxScale;
	vecScale=dblBase.^vecExp;
	intScaleNum = numel(vecScale);
	matMSD = zeros(intN,intScaleNum);
	if boolUseParallel
		parfor intScaleIdx=1:intScaleNum
			dblScale = vecScale(intScaleIdx);
			
			%run through all points
			for intS=1:intN
				%select points within window
				matMSD(intS,intScaleIdx) = getD(dblScale,intS,intN,vecT,vecV);
			end
		end
	else
		for intScaleIdx=1:intScaleNum
			dblScale = vecScale(intScaleIdx);
			
			%run through all points
			for intS=1:intN
				%select points within window
				matMSD(intS,intScaleIdx) = getD(dblScale,intS,intN,vecT,vecV);
			end
		end
	end
	
	%% smoothing
	if intSmoothSd > 0
		vecFilt = normpdf(-2*(intSmoothSd):2*intSmoothSd,0,intSmoothSd)';
		vecFilt = vecFilt./sum(vecFilt);
		%pad array
		matMSD = padarray(matMSD,floor(size(vecFilt)/2),'replicate');
		
		%filter
		matMSD = conv2(matMSD,vecFilt,'valid');
		
		%title
		strTitle = 'Smoothed MSDs';
	else
		%title
		strTitle = 'MSDs';
	end
	%mean
	vecM = mean(matMSD,2);
	
	%weighted average of vecM by inter-spike intervals
	dblMeanM = (1/dblUseMaxDur)*sum(((vecM(1:(end-1)) + vecM(2:end))/2).*diff(vecT));
	
	%rescale to real firing rates
	vecRate = dblMeanRate * ((vecM + 1/dblUseMaxDur)/(dblMeanM + 1/dblUseMaxDur));
		
	%% plot
	if intPlot == 1
		if ~isempty(get(gca,'Children'))
			figure;
		end
		stairs(vecT,vecRate)
		xlabel('Time after event (s)');
		ylabel(strLabelY);
		title(sprintf('Peri Event Plot (PEP)'));
		fixfig
	elseif intPlot > 1
		subplot(2,3,5);
		imagesc(matMSD');
		set(gca,'ytick',[]);
		ylabel(sprintf('Scale (s) (%.1es - %.1es)',vecScale(1),vecScale(end)));
		xlabel('Timestamp index (#)');
		title(strTitle);
		fixfig
		grid off
		
		subplot(2,3,6);
		if numel(vecT) > 10000
			vecSubset = round(linspace(1,numel(vecT),10000));
			plot(vecT(vecSubset),vecRate(vecSubset));
		else
			stairs(vecT,vecRate);
		end
		xlabel('Time after event (s)');
		ylabel(strLabelY);
		title(sprintf('Peri Event Plot (PEP)'));
		fixfig
	end
	
	%% build output
	if nargout > 1
		sMSD = struct;
		sMSD.vecRate = vecRate;
		sMSD.vecT = vecT;
		sMSD.vecM = vecM;
		sMSD.vecScale = vecScale;
		sMSD.matMSD = matMSD;
		sMSD.vecV = vecV;
	end
end
function dblD = getD(dblScale,intS,intN,vecT,vecV)
	%select points within window
	dblT = vecT(intS);
	dblMinEdge = dblT - dblScale/2;
	dblMaxEdge = dblT + dblScale/2;
	intIdxMinT = find(vecT > dblMinEdge,1);
	if isempty(intIdxMinT),intIdxMinT=1;end
	intIdxMaxT = find(vecT > dblMaxEdge,1);
	if isempty(intIdxMaxT),intIdxMaxT=intN;end
	if intIdxMinT == intIdxMaxT && intIdxMinT > 1,intIdxMinT=intIdxMaxT-1;end
	dbl_dT = max([dblScale (vecT(intIdxMaxT) - vecT(intIdxMinT))]);
	dblD = (vecV(intIdxMaxT) - vecV(intIdxMinT))/dbl_dT;
end
