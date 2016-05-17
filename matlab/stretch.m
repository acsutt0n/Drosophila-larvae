% out = stretch(in, newlength, varargin)
%
% Stretch an array so that it's length matches that of another array
%

function out = stretch(in, newlength, varargin)

if nargin < 3
  how = 'interp';
else
  how = varargin{1};
end

out = zeros(newlength, 1);

switch how
  case 'end'
    out(1:length(in)) = in;
    for u = length(in)+1:newlength
      out(u) = in(end);
    end
  
  case 'interp'
    factor = ceil(newlength/length(in));
    if factor > 1 % Else can't really do anything
      out = interp(in, ceil(newlength/length(in)));
      out = out(1:newlength);
    else
      fprintf('Cannot work with a factor of %i', factor);
    end
end



end