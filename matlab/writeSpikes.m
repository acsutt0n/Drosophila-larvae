% writeSpikes(fname, options)
%
% This calculates the moving spikes from the given file.
% It needs access to abfload, GetSpikes and movingSpikes. The output is written
% as a csv file, even though these might be get pretty large. Can't
% think of a better way to do it
%
% INPUT:
%     fname - abf file name
%     options - this is blank for now
%
% OUTPUT - no output is yielded, but a .csv is generated.

function writeSpikes(fname, varargin)

% Inputs
if nargin < 2
  asstruct = true;
else
  asstruct = varargin{1};
end

% File name
splits = strsplit(fname, '.');
splits{2} = '_props.csv';
outfile = strjoin(splits,'');

% First open up the abf file and send its main (only?) trace to GetSpikes
D = abfload(fname);
% Trace should be in location D(:,1), but we check just in case
[~ ,n] = size(D);
use = 1;
for i=1:n
  if mean(D(:,i)) < -10 % Resting potential should be around -30 mV, 
    use = i;            % so this should catch it
  end
end

% Get the spikes
fprintf('Getting the spike struct with a moving window...\n');
spk = movingSpikes(0.1, D(:, 1));

if asstruct
  % Keep it as a struct, but transpose the rows etc
  n_spk = {};
  fields = fieldnames(spk);
  for i=1:length(fields)
    n_spk.(fields{i}) = spk.(fields{i})';
  end
  
  % Save as a struct
  struct2csv(n_spk, outfile);
  fprintf('File %s written.\n', outfile);

else
  % Create a matrix to write to csv
  % Find the max field length
  maxF = 0;
  for i=1:length(fields)
    chk = length(spk.(fields{i}));
    % fprintf('Field: %s %i\n', fields{i}, chk);
    if chk > maxF
      maxF = chk;
    end
  end

  % Create the matrix and populate it
  F = zeros(length(fields), maxF);

  for i=1:length(fields)
    F(i,1:length(spk.(fields{i}))) = spk.(fields{i});
  end
  F(F==0) = nan; % Replace 0 with NaN
  
  % Write to csv file
  csvwrite(outfile, F);
  fprintf('File %s written.\n', outfile);
end


end





