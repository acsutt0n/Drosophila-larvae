% test.m %%%%%%%%%%%% testing stuff

function [peakInds, peakVals] = movingMax(dT, trace, varargin)

% Deal with the inputs first
if nargin < 3
  rise = 2; % mV
else
  rise = varargin{1};
end
if nargin < 4
  run = 5; % ms
else
  run = varargin{2};
end
if nargin < 5
  win = run/dT; % Default is the whole run time for moving window
else
  win = varargin{3};
end

% Calculate the moving windows 
fprintf('Calculating %i moving windows ...', length(trace)-win-1)
diffs = diff(trace);
dVs = [];
for i=1:length(diffs)-win % Scroll through the moving window
  dVs = [dVs; sum(diffs(i:i+win))];
end

% Find some candidate dVs
dVs = [dVs, [1:length(dVs)]];
sorted_dVs = sortrows(dVs, 1);


% Show how these compare to the trace
hold on;
plot(linspace(0,length(trace)/1000*dT, dT), trace)
