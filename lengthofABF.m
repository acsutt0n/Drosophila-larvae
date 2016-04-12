%
% len = lengthofABF(abfname, dT)
%
%  Returns the length, in ms, of the trace.
%
%  INPUT:
%      abfname - path of abf file that you want to check
%      dT - in ms, 0.1 by default
%
%  OUTPUT:
%      len - length of the trace, in ms
%

function len = lengthofABF(abfnames, dT)

% if sequence, load one at a time
if ischar(abfnames)
  abfnames = {abfnames};
end

len = zeros(length(abfnames),1);

% Check input for dT
if nargin < 2
  dT = 0.1;
else
  dT = dT;
end

for p = 1:length(abfnames)
  % Make sure abfload is accessible
  try
    trace = abfload(abfnames{p});
  catch
    addpath(genpath('/home/alex/code/'));
    trace = abfload(abfnames{p});
  end

  % Figure out length of trace
  trace = trace(:,1);
  len(p) = length(trace)*dT;


end









