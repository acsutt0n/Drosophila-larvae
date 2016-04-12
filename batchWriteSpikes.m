% excepts = batchWriteSpikes(filelist)
%
% Apply writeSpikes to all abf files in a list, as a txt file.
%
% INPUT:
%  filelist - a txt file of filepaths, each line is an abf file
%
% OUTPUT:
%  excepts - a list of files that threw exceptions

function excepts = batchWriteSpikes(filelist)

% Get the file list and run them through writeSpikes
flist = importdata(filelist);
excepts = {};

for i=1:length(flist)
  try
    writeSpikes(flist{i})
  catch
    excepts{i} = flist{i};
  end
end



end

