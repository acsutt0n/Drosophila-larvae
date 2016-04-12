% [locs, powers] = spectralSpikes(trace, dT)
%
%  spectralSpikes.m - Spectal spike analysis for bursts. This also smooths
%  but includes the raw in the plot (only returns the smoothed bins).
%  
%  Input:
%      trace - N x 1 list of voltages
%      dT - in ms, defaults to 0.1 (10kHz)
%
%  Output:
%      locs - x-axis, frequency locations
%      powers - y-axis, relative powers of each frequency
%

function [locs, powers] = spectralSpikes(trace, dT)

if nargin < 2
  dT = 0.1;
end

% Get the sampling frequency and other stuff
Fs = 1000/dT; 
N = length(trace);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run the FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%

xdft = fft(trace);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;
dfreq = freq(2)-freq(1); % delta frequency steps
dind = ceil(1/dfreq); % delta index for 1 Hz

% Filter and smooth the data
ppwr = 10*log10(psdx);
dpwr = downsample(ppwr, 100);
mavg = movingAvg(dpwr, 200); % Get a moving avg
mavg = stretch(mavg, length(ppwr), 'interp'); % Resize the moving avg
newp = ppwr - mavg; % Subtract the slow baseline

% Sovitzky-Golay filter: 
sfit = sgolayfilt(ppwr, 3, 111);
sfitavg = movingAvg(sfit, 10);
sfit2 = sfit-mavg;

% Use S-G to find the locations and powers
freq100 = min(find(freq>100));
trough = min(sfit2(1:freq100/10));
trind = find(sfit2==trough); % This is trough

% Find the first peak by slope
slope = diff(sfit2);
movSlope = [];
for s = trind:freq100
  start = max([1, s-floor(dind/2)]);
  fin = min([length(sfit2), s+floor(dind/2)]);
  % The moving slope window is 1 Hz (+/- 0.5 Hz)
  movSlope = [movSlope; mean(sfit2(start:fin))];
end
for m = 

slopeind = find(movSlope==max(movSlope));

% Find the "first" peak
maxpk = max(sfit2(trind:freq100)); % Peak starting at trough
maxind = find(sfit2==maxpk); % This is peak
% If peaks don't agree, make two, else just one
if abs(maxind - slopeind) < dind
  locs = [freq(slopeind); freq(maxind)];
  powers = [sfit2(slopeind); maxpk];
else % Just register one
  locs = [freq(maxind)];
  powers = [maxpk]; 
end


% Fit a decay to this and find another peak > 1 std above the noise
Ps = polyfit(freq(maxind:freq100)', sfit2(maxind:freq100), 2);
pfit = Ps(1)*freq(maxind:freq100).^2 + Ps(2)*freq(maxind:freq100)+Ps(3);
s_pfit = sfit2(maxind:freq100)';
stans = std(s_pfit-pfit); % Std of sfit2-polyfit
pospks = []; % Possible peaks

for s = 1:length(s_pfit)
  if s_pfit(s) - pfit(s) > 3*stans && s > dind % Ensures it's not the first peak
    pospks = [pospks; s];
  end
end

% For each possible peak, at least +/- 1 Hz should also be > 3*stans
newpots = []; % New potential peaks
if length(pospks) > 0
  for pot = 1:length(pospks)
    start = max([1, pot-dind]);
    fin = min([length(s_pfit), pot+dind]);
    positive = s_pfit(start:fin)> 1.5*stans;
    if prod(positive) == 1
      newpots = [newpots; pot];
    end
  end
  % Make sure peaks aren't too close (+/- dind)
end
if length(newpots) > 0
  locs = [locs; freq(newpots(1)+maxind)];
  powers = [powers; sfit2(newpots(1)+maxind)];
  current = 1;
  for pot = 1:length(newpots)
    if newpots(pot) - newpots(current) > dind
      locs = [locs; freq(newpots(pot)+maxind)];
      powers = [powers; sfit2(newpots(pot)+maxind)];
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show everything so far
figure();
hold on;
a1 = plot(freq, newp, 'b'); alpha(0.3);
a2 = plot(freq, sfit2, '--r'); alpha(0.3);
a3 = plot(freq, stretch(sfitavg, length(newp), 'interp')-mavg, '-g' ); alpha(0.3);
a4 = plot(freq(maxind:freq100), pfit, '--k');
set(a2, 'linewidth', 3); set(a3, 'linewidth', 3); set(a4, 'linewidth', 2);

for pk = 1:length(locs)
  plot(locs(pk), powers(pk), 'r*')
end

legend('Mean-subtracted PSD' , 'Sovitzky-Golay smoothed' , 'Moving avg', ...
       'polyFit'); 


end

