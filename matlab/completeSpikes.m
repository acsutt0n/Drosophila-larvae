% spikes = completeSpikes(spikes, trace)
%   Adds a few other fields to the spikes structure. Includes:
%     - rise time
%     - decay time
%     - A.H.P.
%     - Vthresh
% 

function spikes = completeSpikes(spikes, trace)

% Initialize vectors
spikes.riseTime = [];
spikes.decayTime = [];
spikes.AHP = [];
spikes.Vthresh = [];

% Calculate these metrics for each spike...
for i=1:length(spikes.times)
  % Get the temporary waveform
  wave = trace(spikes.n1List(i):spikes.n2List(i));
  
  % Rise time
  







end












