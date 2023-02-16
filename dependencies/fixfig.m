function fixfig(handle,boolMakeActive,dblLineWidth)
	
	%inputs
	if ~exist('handle','var') || isempty(handle)
		handle=gca;
	end
	if ~exist('boolMakeActive','var') || isempty(boolMakeActive)
		boolMakeActive=true;
	end
	if ~exist('dblLineWidth','var') || isempty(dblLineWidth)
		dblLineWidth=2;
	end
	
	%check if axes or figure is supplied
	if isa(handle,'matlab.ui.Figure')
		if isfield(handle.Children(end),'Axes')
			handle = handle.Children(end).Axes;
		else
			error([mfilename ':HandleIsFig'],'Handle is a figure, not Axes');
		end
	end
		
	if boolMakeActive
		axes(handle);
	end
	
	%remove box
	set(handle, 'box','off');
	drawnow;
	dblFontSize=14; %change the figure font size
	try 
		xlabel(get(get(handle,'xlabel'), 'String'),'FontSize',dblFontSize); %set x-label and change font size
		ylabel(get(get(handle,'ylabel'), 'String'),'FontSize',dblFontSize);%set y-label and change font size
	catch
		%not a cartesian plot
	end
	title(get(get(handle,'title'),'string'),'FontSize',14);
	set(handle,'FontSize',dblFontSize,'Linewidth',dblLineWidth); %set grid line width and change font size of x/y ticks
	set(handle,'TickDir', 'out');
	if ~strcmp(get(get(handle,'Children'),'Type'),'image')
		vecChildren = get(handle,'Children');
		for intChild=1:numel(vecChildren)
			if isprop(vecChildren(intChild),'Linewidth')
				set(vecChildren(intChild),'Linewidth',dblLineWidth);end %change default linewidth to 2
		end
	end
end

