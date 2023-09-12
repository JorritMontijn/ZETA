function fixfig(handle,boolMakeActive,dblLineWidth,dblFontSize)
	
	%flush
	drawnow;
	pause(1e-3);
	
	%inputs
	if ~exist('handle','var') || isempty(handle)
		handle=gcf;
	end
	if ~exist('boolMakeActive','var') || isempty(boolMakeActive)
		boolMakeActive=true;
	end
	if ~exist('dblLineWidth','var') || isempty(dblLineWidth)
		dblLineWidth=2;
	end
	if ~exist('dblFontSize','var') || isempty(dblFontSize)
		dblFontSize=14; %change the figure font size
	end
	
	%make active
	if boolMakeActive && isaxes(handle)
		axes(handle);
	elseif boolMakeActive && isa(handle,'matlab.ui.Figure')
		figure(handle);
	end
	
	%check if axes or figure is supplied
	if isa(handle,'matlab.ui.Figure')
		for intChild=1:numel(handle.Children)
			if isfield(handle.Children(intChild),'Axes')
				fixfig(handle.Children(intChild).Axes,boolMakeActive,dblLineWidth,dblFontSize);
			elseif isfield(handle.Children(intChild),'Colorbar')
				fixfig(handle.Children(intChild).Colorbar,boolMakeActive,dblLineWidth,dblFontSize);
			elseif isaxes(handle.Children(intChild))
				fixfig(handle.Children(intChild),boolMakeActive,dblLineWidth,dblFontSize);
			elseif isa(handle.Children(intChild),'matlab.graphics.illustration.ColorBar')
				fixfig(handle.Children(intChild),boolMakeActive,dblLineWidth,dblFontSize);
			end
		end
		return;
	end
	
	%remove box
	if isaxes(handle)
		set(handle, 'box','off');
	end
	
	drawnow;
	try
		set(get(handle,'xlabel'),'FontSize',dblFontSize); %set x-label and change font size
		set(get(handle,'ylabel'),'FontSize',dblFontSize);%set y-label and change font size
	catch
		%not a cartesian plot
	end
	set(handle,'FontSize',dblFontSize,'Linewidth',dblLineWidth); %set grid line width and change font size of x/y ticks
	set(handle,'TickDir', 'out');
	if ~strcmp(get(get(handle,'Children'),'Type'),'image')
		vecChildren = get(handle,'Children');
		for intChild=1:numel(vecChildren)
			if isprop(vecChildren(intChild),'Linewidth')
				set(vecChildren(intChild),'Linewidth',dblLineWidth);end %change default linewidth to 2
		end
	end
	drawnow;
end

