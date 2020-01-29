function c = redbluepurple(m)
	%REDBLUEPURPLE
	
	if nargin < 1, m = size(get(gcf,'colormap'),1); end
	
	% From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
	r = linspace(0,1,m)';
	g = 0*r;
	b = flipud(r);
	
	c = [r g b];
	
