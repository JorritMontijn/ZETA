function [handleFill,handleLine] = errorfill(vecX,vecY,vecErr1,varargin)
	%errorfill Plot mean +- shaded error region
	%   Syntax: [handleFill,handleLine] = errorfill(vecX,vecY,vecErr,[vecErr2],[vecColor])
	%
	%	input:
	%	- vecX; vector with x-values
	%	- vecY; vector with y-values
	%	- vecErr; vector with error-values for shaded area
	%	- [vecErr2]; optional, if supplied the above vector is used as
	%				upper-limits, and this vector is used as lower-limits
	%	- vecColor; optional 3-element RGB colour vector (default is blue)
	%	
	%	Version history:
	%	1.0 - August 1 2013
	%	Created by Jorrit Montijn
	%	2.0 - Feb 6 2019
	%	Updated to use alpha-mapping for transparency of shaded area
	
	
	%% check inputs
	if ~isempty(varargin) && length(varargin{1}) == length(vecErr1)
		intArgOffset = 1;
		vecErr2 = varargin{1};
	else
		intArgOffset = 0;
		vecErr2 = vecErr1;
	end
	%% prep inputs
	%switch orientation
	if isvector(vecX)
		vecX = vecX(:);
	end
	if isvector(vecY)
		vecY = vecY(:);
		vecErr1 = vecErr1(:);
		vecErr2 = vecErr2(:);
	end
	
	%% check color
	if nargin >= 4+intArgOffset && length(varargin{1+intArgOffset}) == 3
		vecColor = varargin{1+intArgOffset};
	else
		vecColor = lines(size(vecX,2));
	end
	
	%get selections
	intX = length(vecX);
	vecWindowInv = intX:-1:1;
	vecXinv = vecX(vecWindowInv,:);
	
	%plot lines first to allow legend to show lines only
	hold on
	for intCurve = 1:size(vecX,2)
		handleLine = plot(vecX(:,intCurve),vecY(:,intCurve),'-','LineWidth',2,'Color',vecColor(intCurve,:));
	end
	for intCurve = 1:size(vecX,2)
		handleFill = patch([vecX(:,intCurve);vecXinv(:,intCurve)],[vecY(:,intCurve)+vecErr1(:,intCurve,1);vecY(vecWindowInv,intCurve)-vecErr2(vecWindowInv,intCurve,end)],vecColor(intCurve,:),'EdgeColor','none');
		alpha(handleFill,.5);
	end
	hold off
	drawnow;
end

