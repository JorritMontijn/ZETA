function output = roundi(number,digits,argument)
	%roundi rounds a number to a certain amount of digits
	%syntax:
	% output = roundi(number,digits,argument)
	%
	%input variables:
	% number: the input number
	% digits: required amount of digits
	% OPTIONAL: argument: 'floor', 'ceil' or 'fix'; this makes the function use not the round function, but the function specified by argument
	% By Jorrit Montijn [21-10 2009]
	
	temp = number*(10^(digits));
	if exist('argument','var')
		if strcmp(argument,'floor')
			rounded = floor(temp);
		elseif strcmp(argument,'ceil')
			rounded = ceil(temp);
		elseif strcmp(argument,'fix')
			rounded = fix(temp);
		end
	else
		rounded = round(temp);
	end
	output = rounded/(10^(digits));
end
