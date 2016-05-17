%
%  bursttms = idBursts(paths, dT)
%
%

function [bursttms, failed] = idBursts(paths, dT, varargin)

if nargin < 2
  dT = 0.0001; % in s!
else
  dT = dT; % in s!
end

bursttms = {};
failed = {};

for fn = 1:length(paths)
  % Load the trace
  trace = abfload(paths{fn});
  trace = trace(:,1);
  
  if nargin > 2
    gg = varargin{1};
    keepstart = floor( gg(fn,1) / dT);
    keepstop = floor( gg(fn,2) / dT);
  else
    keepstart = 1;
    keepstop = length(trace);
  end
  
  
  fprintf('start: %i, stop %i\n', keepstart, keepstop);
  
  if keepstop < 1 || keepstop > length(trace)
    keepstop = length(trace);
    fprintf('Tried to start at %i and STOP at %i\n', keepstart, keepstop);
  end
  if keepstart < 1 || keepstart > length(trace)
    fprintf('Tried to START at %i and stop at %i\n', keepstart, keepstop);
    keepstart = 1;
  end
  trace = trace(keepstart: keepstop); % Cut the trace if needed
  
  % Display the trace for manual burst targeting
  tms = dT:dT:length(trace)*dT;
  split = 1;
  if length(tms)*dT > 60 % 1 min
    split = 2;
  elseif length(tms)*dT > 120 % 2 min
    split = 3;
  else
    split = 4;
  end
  
  for u = 1:split
    sta = (length(trace)/split)*(u-1) + 1;
    sto = (length(trace)/split)*(u);
    subplot(split,1,u), plot(tms(sta:sto), trace(sta:sto));
  end
  
  title('Identify BURSTS! Press Enter to stop. One RIGHT CLICK indicates failed burst');
  [x,~, button] = ginput(1000);
  
  
  % Identify and record the burst windows
  try
    tfil = strsplit(paths{fn}, '/');
    tfil = tfil(end);
    tfil = strsplit(char(tfil), '.');
    tfil = tfil(1);
  catch
    tfil = paths(fn);
  end
  
  tfil = strcat('id_', tfil);
  
  % fprintf('%s', tfil);
  durs = []; fails = [];
  cnt = 1;
  for i = 1:length(x)
    if button(i) == 1
      if mod(cnt,2) == 0
        durs(cnt/2,2) = x(i);
      else
        durs = [durs; x(i), 0];
      end
      cnt = cnt + 1;
    elseif button(i) == 3 % Right mouse click indicates a failed burst
      fails = [fails; x(i)];
    end
  end
  
  % Write into the structs
  disp(char(tfil));
  bursttms = setfield(bursttms, char(tfil), durs);
  failed = setfield(failed, char(tfil), fails);
end


  












end




