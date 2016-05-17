% spikes = movingSpikes(dT, trace, varargin)
%   This runs GetSpikes over multiple windows and compiles the results.
%   If there is a descrepency between spike voltages, for example if the 
%   recording deteriorates, this will still be able to identify spikes
%   when GetSpikes alone fails.
%
% INPUT
%   dT: time step of recording; if > 1 then it's inverse is taken (in ms).
%   trace: a 1-D vector of voltage recordings.
%   numWindows: number of windows over which GetSpikes is called;
%               defaults to 1 per second.
%
% OUTPUT
%   spikes: a struct of compiled spike times, etc; as per GetSpikes:
%    -spike.times is a plain list of spike times (in ms)
%    -spike.height is a plain list of spike heights (in mV)
%    -spike.width is a plain list of spike width (in ms)
%    -spike.freq is overall spiking frequency (in Hz)
%    -spike.intervals is a list of interspike intervals (in ms)
%    -spike.frequencies is a list of instantaneous frequencies (in Hz)
%    Shape information structures (should be self-descriptive)
%    -spike.maxV, spike.maxDeriv, spike.minDeriv, spike.preMinV,
%     spike.postMinV, spike.preMaxCurve, spike.postMaxCurve
%           Each contains a list of times/voltage points, and if relevant
%           another quantity (such as K for curvatures)
%


function spikes = movingSpikes(dT, trace, varargin)

if nargin < 3
  numWindows = length(trace) * dT / 1000; % 1 window per second
else
  numWindows = varargin{1}; % Throws an error if something funky happens
end

% Set the defaults
spikes = {};
spikes.times = [];
spikes.n1List = []; % These are pulled from n1List:n2List
spikes.n2List = []; 
spikes.maxV = [];
spikes.maxVtms = [];
spikes.maxVinds = []; 
spikes.maxDerivV = []; % Max deriv
spikes.maxDerivdV = [];
spikes.maxDerivtms = [];
spikes.maxDerivinds = [];
spikes.minDerivV = []; % Min deriv
spikes.minDerivdV = [];
spikes.minDerivtms = [];
spikes.minDerivinds = [];
spikes.preMinV = []; % pre-spike min voltage
spikes.preMintms = [];
spikes.preMininds = [];
spikes.postMinV = []; % post-spike min voltage
spikes.postMintms = [];
spikes.postMininds = [];
spikes.preMaxCurveV = []; % pre-max curve
spikes.preMaxCurveK = [];
spikes.preMaxCurvetms = [];
spikes.preMaxCurveinds = [];
spikes.postMaxCurveV = []; % post-max curve
spikes.postMaxCurveK = [];
spikes.postMaxCurvetms = [];
spikes.postMaxCurveinds = [];
spikes.height = []; % More data
spikes.repolarizationV = [];
spikes.intervals = []; % ms between spikes
spikes.frequencies = []; % Instantaneous frequencies
spikes.freq = []; % Avg freqs (one per window)
% spikes.winMid = []; % The middle of the current moving window


for s=1:numWindows-1
  % Figure out which interval to analyze
  interval = length(trace) / numWindows;
  if s ~= numWindows % As long as it's not the last one, treat it normally
    start_int = floor(s*interval);
    stop_int = floor((s+1)*interval); % These need to be integers
  else
    start_int = floor(s*interval);
    stop_int = length(trace);
  end
  
  if mod(s,50) == 0
    fprintf('Window %i of %i ....\n', s, numWindows)
  end
  
  % Collect spk features for this sub-trace
  try
    spk = GetSpikes(dT, trace(start_int:stop_int));
  catch
    length(trace)
    start_int
    spk = GetSpikes(dT, trace(start_int:length(trace)));
  end
  % Populate the data fields
  spikes = populateSpikeStruct(spikes, spk, dT, start_int);
  
  
end

% Now have all the spikes data.
spikes;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populate the structure spikes with new data from GetSpikes
%  dT and indOffset are used to offset spike times and indices.
%  indOffset: start_int above.

function spikes = populateSpikeStruct(spikes, spk, dT, indOffset)
  
spikes.times = [spikes.times, spk.times+dT*indOffset];
spikes.n1List = [spikes.n1List, spk.n1List+indOffset]; % These are used to re-create
spikes.n2List = [spikes.n2List, spk.n2List+indOffset]; % the individual waveforms
spikes.maxV = [spikes.maxV, spk.maxV.v];
spikes.maxVtms = [spikes.maxVtms, spk.maxV.t+dT*indOffset];
spikes.maxVinds = [spikes.maxVinds, spk.maxV.ind+indOffset]; 
spikes.maxDerivV = [spikes.maxDerivV, spk.maxDeriv.v]; % Max deriv
spikes.maxDerivdV = [spikes.maxDerivdV, spk.maxDeriv.dV];
spikes.maxDerivtms = [spikes.maxDerivtms, spk.maxDeriv.t+dT*indOffset];
spikes.maxDerivinds = [spikes.maxDerivinds, spk.maxDeriv.ind+indOffset];
spikes.minDerivV = [spikes.minDerivV, spk.minDeriv.v]; % Min deriv
spikes.minDerivdV = [spikes.minDerivdV, spk.minDeriv.dV];
spikes.minDerivtms = [spikes.minDerivtms, spk.minDeriv.t+dT*indOffset];
spikes.minDerivinds = [spikes.minDerivinds, spk.minDeriv.ind+indOffset];
spikes.preMinV = [spikes.preMinV, spk.preMinV.v]; % pre-spike min voltage
spikes.preMintms = [spikes.preMintms, spk.preMinV.t+dT*indOffset];
spikes.preMininds = [spikes.preMininds, spk.preMinV.ind+indOffset];
spikes.postMinV = [spikes.preMinV, spk.preMinV.v]; % post-spike min voltage
spikes.postMintms = [spikes.preMintms, spk.preMinV.t+dT*indOffset];
spikes.postMininds = [spikes.preMininds, spk.preMinV.ind+indOffset];
spikes.preMaxCurveV = [spikes.preMaxCurveV, spk.preMaxCurve.v]; % pre-max curve
spikes.preMaxCurveK = [spikes.preMaxCurveK, spk.preMaxCurve.K];
spikes.preMaxCurvetms = [spikes.preMaxCurvetms, spk.preMaxCurve.t+dT*indOffset];
spikes.preMaxCurveinds = [spikes.preMaxCurveinds, spk.preMaxCurve.ind+indOffset];
spikes.postMaxCurveV = [spikes.postMaxCurveV, spk.postMaxCurve.v]; % post-max curve
spikes.postMaxCurveK = [spikes.postMaxCurveK, spk.postMaxCurve.K];
spikes.postMaxCurvetms = [spikes.postMaxCurvetms, spk.postMaxCurve.t+dT*indOffset];
spikes.postMaxCurveinds = [spikes.postMaxCurveinds, spk.postMaxCurve.ind+indOffset];
spikes.height = [spikes.height, spk.height]; % More data
spikes.repolarizationV = [spikes.repolarizationV, spk.repolarizationPotential];
spikes.intervals = [spikes.intervals, spk.intervals]; % ms between spikes
spikes.frequencies = [spikes.frequencies, spk.intervals]; % Instantaneous frequencies
spikes.freq = [spikes.freq, spk.freq]; % Avg freqs (one per window)
  
  
end