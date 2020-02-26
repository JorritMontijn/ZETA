function [dblOnset,dblValue] = getOnset(vecData,vecT,dblPeakT,vecRestrictRange)
	%getOnset Returns peak onset. Syntax:
	%    [dblOnset,dblValue] = getOnset(vecData,vecT,dblPeakT,vecRestrictRange)
	%
	%Required input:
	%	- vecData [N x 1]: values
	%
	%Optional inputs:
	%	- vecT [N x 1]: timestamps corresponding to vecData (default: [1:N])
	%	- dblPeakT (float): timestamp corresponding to peak
	%	- vecRestrictRange [2 x 1]: restrict peak to lie within vecRestrictRange(1) and vecRestrictRange(end)
	%
	%Outputs:
	%	- dblOnset: time of peak onset (first crossing half-height of peak)
	%	- dblValue: value at peak onset
	%
	%Version history:
	%1.0 - February 26 2020
	%	Created by Jorrit Montijn

	%check inputs
	if ~exist('vecT','var')
		vecT = 1:numel(vecData);
	end
	if ~exist('vecRestrictRange','var')
		vecRestrictRange = [min(vecT) min(vecT)+range(vecT)/2];
	end
	if ~exist('dblPeakT','var')
		[dummy,dblPeakT] = getPeak(vecData,vecT,vecRestrictRange);
	end
	%remove time points outside restricted range
	indRemove = vecT < vecRestrictRange(1) | vecT > vecRestrictRange(end);
	vecT = vecT(~indRemove);
	vecData = vecData(~indRemove);
	
	%calculate first timepoint crossing half-height of peak 
	[dummy,intPeakIdx]=min(abs(vecT-dblPeakT));
	dblPeakVal = vecData(intPeakIdx);
	dblBaseVal = vecData(1);
	dblThresh = (dblPeakVal - dblBaseVal)/2 + dblBaseVal;
	intOnsetIdx = find(vecData > dblThresh,1,'first');
	dblOnset = vecT(intOnsetIdx);
	dblValue = vecData(intOnsetIdx);
	
	%check if empty
	if isempty(intOnsetIdx)
		dblOnset = nan;
		dblValue = nan;
	end
end
