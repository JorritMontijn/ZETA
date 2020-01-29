function [vecMSD,sMSD] = getMultiScaleDeriv(vecT,vecV,intSmoothSd,dblMinScale,dblBase,intPlot)
	%getMultiScaleDeriv Returns multi-scale derivative. Syntax:
	%   [vecMSD,sMSD] = getMultiScaleDeriv(vecX,vecY,intSmoothSd,dblMinScale,dblBase)
	%Required input:
	%	- vecT [N x 1]: timestamps (e.g., spike times)
	%	- vecV [N x 1]: values (e.g., z-scores)
	%
	%Optional inputs:
	%	- intSmoothSd: Gaussian SD of smoothing kernel (in # of bins) [default: 3]
	%	- dblMinScale: minimum derivative scale in seconds [default: 1/1000]
	%	- dblBase: critical value for locally dynamic derivative [default: 4]
	%	- intPlot: integer, plotting switch (0=none, 1=plot)
	%
	%Outputs:
	%	- vecMSD; Multi-scale derivative
	%	- sMSD; structure with fields:
	%		- vecMSD; Multi-scale derivative
	%		- vecScale; timescales used to calculate derivatives
	%		- matSmoothMSD; smoothed multi-scale derivatives matrix
	%		- matMSD; raw multi-scale derivatives matrix
	%
	%Version history:
	%1.0 - January 24 2019
	%	Created by Jorrit Montijn - split from previous getMultiScaleDeriv.
	
	%% set default values
	if ~exist('intSmoothSd','var') || isempty(intSmoothSd)
		intSmoothSd = 5;
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
	
	%% prepare data
	dblMaxScale = log(max(vecT)) / log(dblBase);
	intN = numel(vecT);
	
	%% get multi-scale derivative
	vecExp = dblMinScale:dblMaxScale;
	vecScale=dblBase.^vecExp;
	intScaleNum = numel(vecScale);
	matDeriv = zeros(intN,intScaleNum);
	try %try parallel
		parfor intScaleIdx=1:intScaleNum
			dblScale = vecScale(intScaleIdx);
			
			%run through all points
			for intS=1:intN
				%select points within window
				dblT = vecT(intS);
				dblMinEdge = dblT - dblScale/2;
				dblMaxEdge = dblT + dblScale/2;
				intIdxMinT = find(vecT > dblMinEdge,1);
				if isempty(intIdxMinT),intIdxMinT=1;end
				intIdxMaxT = find(vecT > dblMaxEdge,1);
				if isempty(intIdxMaxT),intIdxMaxT=intN;end
				matDeriv(intS,intScaleIdx) = (vecV(intIdxMaxT) - vecV(intIdxMinT))/dblScale;
			end
		end
	catch %otherwise try normal loop
		for intScaleIdx=1:intScaleNum
			dblScale = vecScale(intScaleIdx);
			
			%run through all points
			for intS=1:intN
				%select points within window
				dblT = vecT(intS);
				dblMinEdge = dblT - dblScale/2;
				dblMaxEdge = dblT + dblScale/2;
				intIdxMinT = find(vecT > dblMinEdge,1);
				if isempty(intIdxMinT),intIdxMinT=1;end
				intIdxMaxT = find(vecT > dblMaxEdge,1);
				if isempty(intIdxMaxT),intIdxMaxT=intN;end
				matDeriv(intS,intScaleIdx) = (vecV(intIdxMaxT) - vecV(intIdxMinT))/dblScale;
			end
		end
	end
	
	%% smoothing
	if intSmoothSd > 0
		vecFilt = normpdf(-2*(intSmoothSd):2*intSmoothSd,0,intSmoothSd)';
		vecFilt = vecFilt./sum(vecFilt);
		matSmoothDeriv = conv2(matDeriv,vecFilt,'same');
	else
		matSmoothDeriv = matDeriv;
	end
	%mean
	vecMSD = mean(matSmoothDeriv,2);
		
	%% plot
	if intPlot > 0
		subplot(2,3,5);
		imagesc(matSmoothDeriv');
		vecTicksY = get(gca,'ytick');
		cellTicksY = cellfun(@sprintf,cellfill('%.1e',size(vecTicksY)),vec2cell(vecScale(vecTicksY)),'UniformOutput',false);
		set(gca,'yticklabel',cellTicksY);
		ytickangle(gca,75)
		ylabel(sprintf('Scale (s) (%.1es - %.1es)',vecScale(1),vecScale(end)));
		xlabel('Timestamp index (#)');
		title('Smoothed MSDs');
		fixfig
		grid off
		
		subplot(2,3,6);
		plot(vecT,vecMSD)
		xlabel('Time (s)');
		ylabel('Time-locked activation');
		title(sprintf('Peri Event Plot (PEP)'));
		fixfig
	end
	
	%% build output
	if nargout > 1
		sMSD = struct;
		sMSD.vecMSD = vecMSD;
		sMSD.vecScale = vecScale;
		sMSD.matSmoothMSD = matSmoothDeriv;
		sMSD.matMSD = matDeriv;
		sMSD.vecT = vecT;
		sMSD.vecV = vecV;
	end
end

