% avg = movingAvg(x, winlength, varargin)
%
%  Computes a moving average, agnostic to step size (dT).
%
%  Input:
%       x - a trace of data
%       winlength - the window length, if = 1 no filtering occurs.
%       toend - the return trace is missing winlength points (due to
%               taper), this allows for two options:
%              'stretch' - stretch the beginning and end 10% to fill the
%                          space
%              'fill'    - simply repeats the first and last to the
%                          beginning/end 
%
%  Output:
%       avg - a trace that has been averaged with a moving window. If
%             'toend' is left blank, avg = length(x) - winlength


function avg = movingAvg(x, winlength, varargin)


if nargin > 2
  toend = varargin{1};
else
  toend = false;
end

% Get the moving window
half = floor(winlength/2);
avg = [];
for u = half+1:length(x)-(half+1)
  avg = [avg; mean(x(u-half:u+half))];
end

% Deal with the tapers
switch toend
  case 'fill'
    newavg = zeros(length(x),1);
    newavg(1:half) = avg(1);
    newavg(half+1:end-half) = avg;
    newavg(end-half:end) = avg(end);
    avg = newavg;
  
  case 'stretch'
    % Stretch is will 10% of beginning/end
    newavg = zeros(length(x),1);
    avg10 = floor(length(avg)*.1);
    mod10 = floor(avg10/half); % Frequency of stretching
    cnt = 1;  %%%%%% First 10%
    for h = 1:half+avg10 
      if mod(h, mod10) == 0 % Repeat the last one every mod10-th value
        newavg(h) = newavg(h-1);
      else % Take the cnt-th value
        newavg(h) = avg(cnt);
        cnt = cnt + 1;
      end
    end
    
    %%%%%% Middle 80%
    newavg(avg10+half+1:length(x)-(avg10+half+1)) = avg(avg10:end-avg10-1);
    
    cnt = 1; %%%%% Last 10%
    for h = length(x)-(half+avg10):length(x) 
      if mod(h, mod10) == 0
        newavg(h) = newavg(h-1);
      else
        newavg(h) = avg(cnt);
        cnt = cnt + 1;
      end
    end
    avg = newavg;
    
  case false % No taper adjustment requested
    fprintf('No _toend_ parameter was given; moving win will be');
    fprintf(' %i (%i shorter than original -- %i\n', length(avg), winlength, length(x));
end














end