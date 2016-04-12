% Get structure of electrophysiology data from .csv file and do early stage
% of analysis on the ephys data
function ephysData = GetEphysData(varargin)
  % define passable options and default values
  parser = inputParser();
  parser.addParameter('masterFile', '')
  parser.addParameter('ephysPath', '/Users/Maria/Desktop/Ephys')
  parser.addParameter('maxChange_Cm', inf)
  parser.addParameter('maxChange_Rm', inf)
  parser.addParameter('maxChange_Ra', inf)
  parser.addParameter('maxChange_Tau', inf)
  parser.addParameter('maxChange_Hold', inf)
  parser.addParameter('keepDataWithNoStats', true)
  
  % get the options
  parser.parse(varargin{:})
  options = parser.Results;
  
  if isempty( options.masterFile )
    ephysPath = options.ephysPath;
    if ephysPath(end) == filesep()
      ephysPath = ephysPath(end-1);
    end
    %ind = find(ephysPath == filesep(), 1, 'last');
    ind = strfind( ephysPath, filesep() ); ind = ind(end);
    pathSpec = ephysPath(1:ind-1);
  elseif exist( options.masterFile, 'dir' )
    pathSpec = options.masterFile;
  else
    pathSpec = '.';
  end
  
  if exist( options.masterFile, 'file' ) ~= 2
    [fileName, pathName] = uigetfile( fullfile(pathSpec, '*.csv'), ...
                                      'Select master file' );
    if ~ischar( fileName )
      ephysData = [];
      return
    end
    options.masterFile = fullfile( pathName, fileName );
  end
  
  % get the information contained in the master file
  masterInfo = getMasterInfo( options );
  
  % convert that information into a structure with information about
  % experiment files
  ephysData = raw2struct( masterInfo );
  
  % get the appropriate analysis for each file
  ephysData = getAnalysis(ephysData, options);
  
  % select experiments where stats don't change too much
  %ephysData = selectStableStats(ephysData, options);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the information about experiments from the spreadsheet
function masterInfo = getMasterInfo(options)
  [~, ~, fileType] = fileparts( options.masterFile );
  if any( strcmp( fileType, {'.xls', '.xlsx', '.ods'} ))
    [~, sheets] = xlsfinfo( options.masterFile );
    [~, ~, sheet] = xlsread( options.masterFile, sheets{1} );
  elseif strcmp( fileType, '.csv' )
    fIn = fopen( options.masterFile, 'r' );
    sheet = {};
    row = 1;
    textLine = fgetl( fIn );
    while ischar( textLine )
      newCells = strsplit( textLine, ',', 'CollapseDelimiters', false );
      sheet(row,1:numel( newCells )) ...
        = cellfun( @(s) strtrim(s), newCells, 'UniformOutput', false); %#ok<AGROW>
      textLine = fgetl( fIn );
      row = row + 1;
    end
    fclose( fIn );
  end
  
  masterInfo = sheet;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the more human-understandable spreadsheet into a more
% machine-understandable structure
function ephysData = raw2struct( sheet )
  % first find the row of data headers.
  % It starts with 'cell ID' in column 1
  cellIdInd = find( strcmpi( sheet, 'cell ID' ), 1 );
  [headerRow, ~] = ind2sub( size( sheet ), cellIdInd );
  %headerRow = find( cellfun(@(c) strcmpi(c, 'cell ID'), sheet(:,1)) );
  headers = sheet(headerRow, :);

  % Figure out the experiment structure by looking at the two rows above
  % the headers
  phases = getExperimentPhases( sheet(headerRow-2,:), headers )
  protocols = getExperimentProtocols( sheet(headerRow-1,:), headers );
  
  % next figure out which columns to use. Nothing to the right of a cell
  % that starts with 'Notes', and nothing that's empty
  notesCol = find( cellfun( @(c) any( strfind( lower( c ), 'notes' ) ), ...
                            headers ) );
  if isempty( notesCol )
    lastCol = numel( headers );
  else
    lastCol = notesCol - 1;
  end
  useCols = arrayfun( @(n) ~isempty( headers{n} ), 1:lastCol );
  
  % only use the information from those columns
  headers = headers(useCols);
  phases = phases(useCols);
  protocols = protocols(useCols);
  
  % read the data into a structure
  ephysData = [];
  numRows = size( sheet, 1 );
  for rowNum = headerRow+1:numRows
    row = sheet(rowNum,useCols);
    
    cellInfo = struct();
    for col = 1:numel( row )
      header = headers{col}; phase = phases{col};
      protocol = protocols{col};
      sHeader = getSafeName(header); sPhase = getSafeName(phase);
      sProtocol = getSafeName(protocol);
      data = row{col};
      if isempty( phase )
        assert( isempty( protocol ), ...
          ['Don''t understand protocol ', protocol, 'with no phase'] )
        cellInfo.(sHeader) = data;
      elseif isempty( protocol )
        cellInfo.(sPhase).(sHeader).fileName = data;
      else
        cellInfo.(sPhase).(sProtocol).(sHeader).fileName = data;
      end
    end
    
    if ~isfield( cellInfo, 'folder' )
      if ispc
        cellInfo.folder = cellInfo.pc_path;
      elseif ismac
        cellInfo.folder = cellInfo.mac_path;
      else
        error( 'Don''t have linux path' )
      end
    end
    ephysData = [ephysData, cellInfo]; %#ok<AGROW>
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure out what conditions (phases) the experiment is divided into, and
% how individual columns fit into the phases
function phases = getExperimentPhases(row, headers)
  phases = cell( 1, numel( headers ) );
  lastNonEmpty = '';
  for n = 1:numel( headers )
    if n <= numel( row ) &&  ~isempty( row{n} )
      lastNonEmpty = row{n};
    end
    phases{n} = lastNonEmpty;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure out what experimental protocols were used, and how individual
% columns fit into those protocols
function protocols = getExperimentProtocols(row, headers)
  protocols = cell(1, length(headers));
  lastNonEmpty = '';

  stats = 'stat'; statsLen = numel( stats );
  for n = 1:numel( headers )
    if n <= numel( row )
      if ~isempty( row{n} )
        lastNonEmpty = row{n};
      elseif numel( headers{n} ) >= statsLen ...
          && strcmpi( headers{n}(1:statsLen), stats )
        lastNonEmpty = '';
      end
    end
    protocols{n} = lastNonEmpty;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get rid of characters matlab doesn't like in variable names
function fieldName = getSafeName(sheetName)
  fieldName = strrep( strrep( sheetName, ' ', '_' ), '-', '' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the appropriate analysis for each experimental file
function ephysData = getAnalysis(ephysFileData, options)
  ephysData = [];
  for n = 1:numel( ephysFileData )
    ephysData = [ephysData, getCellAnalysis( ephysFileData(n), options )]; %#ok<AGROW>
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the appropriate analysis for an experimental file, or traverse the
% files of a given cell and get the appropriate analysis
function cellInfo = getCellAnalysis(cellFileInfo, options, folderName, ...
                                    structName)
if nargin < 3 || isempty(folderName)
  if isfield(cellFileInfo, 'folder')
    % replace Ephys with the correct full path
    folderName = strrep( cellFileInfo.folder, 'Ephys', options.ephysPath );
    % make sure that filesep is consistant
    folderName = strrep( strrep( folderName, '\', filesep ), ...
                         '/', filesep );
  else
    folderName = '';
  end
end

if isfield(cellFileInfo, 'fileName')
  % this is the lowest level of the data structure, get the data + analysis
  fileName = cellFileInfo.fileName;
  if isempty( fileName )
    fileType = '';
  else
    [~,~,fileType] = fileparts(fileName);
  end
  switch fileType
    case '.sta'
      cellInfo = getAtfAnalysis(folderName, fileName, structName, options);
    case '.abf'
      cellInfo = getAbfAnalysis(folderName, fileName, structName, options);
    case ''
      % the file was empty
      cellInfo = struct();
    otherwise
      error('Error analyzing %s - Don''t know how to open files of type %s', fileName, fileType)
  end
else
  % need to traverse further in the structure
  cellInfo = struct();
  for fieldName = fieldnames(cellFileInfo)'
    fieldName = fieldName{1}; %#ok<FXSET>
    if isstruct(cellFileInfo.(fieldName))
      cellInfo.(fieldName) = getCellAnalysis(cellFileInfo.(fieldName), ...
                                           options, folderName, fieldName);
    else
      cellInfo.(fieldName) = cellFileInfo.(fieldName);
    end
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analysis = getAtfAnalysis(folderName, fileName, structName, options) %#ok<INUSD>
% load data from .atf file and return it
if ~exist(folderName, 'dir')
  folders = strsplit(folderName, filesep);
  if isempty(folders{end})
    folderName = folders{end-1};
  else
    folderName = folders{end};
  end
end

if exist(folderName, 'dir') && ~isempty(fileName)
  data = LoadAtf( fullfile(folderName, fileName) );
  analysis = data;
else
  analysis = struct();
end
analysis.fileName = fileName;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analysis = getAbfAnalysis(folderName, fileName, structName, options) %#ok<INUSD>
% load data from .abf file and perform appropriate analysis
if ~exist(folderName, 'dir')
  folders = strsplit(folderName, filesep);
  if isempty(folders{end})
    folderName = folders{end-1};
  else
    folderName = folders{end};
  end
end
if exist(folderName, 'dir') && ~isempty(fileName)
  try
    fullName = fullfile( folderName, fileName );
    fprintf('Analyzing %s\n', fullName)
    data = LoadAbf( fullName );
    if strcmpi(structName, 'ISteps') || ...
        (numel( structName ) >= 2 && strcmpi( structName(1:2), 'FI' ) )
      close all
      analysis = GetFICurve(data, 'plotSubject', fullName, 'debugPlots', true);
      analysis.data = data;
    elseif strcmpi(structName, 'Ramp')
      analysis = AnalyzeCurrentRamp(data);
    else
      error('Don''t know how to analyze %s protocol', structName)
    end
  catch err
    causeID = 'GetEphysData:getAbfAnalysis';
    causeMessage = sprintf('Error processing %s', ...
                           fullfile(folderName, fileName));
    causeException = MException( causeID, causeMessage );
    rethrow( addCause( err, causeException ) )
  end
else
  analysis = struct();
end
analysis.fileName = fileName;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ephysData = selectStableStats(ephysData, options)
% select experiments where stats don't change too much
for n = 1:length(ephysData)
  cell = ephysData(n);
  for field = fieldnames(cell)'
    field = field{1}; %#ok<FXSET>
    if isstruct(cell.(field))
      % this should be a struct containing stats and other files
      ephysData(n).(field) = selectStableCellStats(cell.(field), options);
    end
  end
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = selectStableCellStats(results, options)
fieldNames = fieldnames(results);
isStat = cellfun( @(s) any(strfind(s, 'stats')), fieldNames);
assert(sum(isStat) == 2, ...
       'Don''t know how to interpret conditions without two stats files')
statsA = results.(fieldNames{find(isStat, 1)});
statsB = results.(fieldNames{find(isStat, 1, 'last')});
if isfield(statsA, 'data')
  statsA = statsA.data;
elseif options.keepDataWithNoStats
  return % statsA is missing, so keep origResults
else
  results = struct(); % statsA is missing, so call results unstable
  return
end
if isfield(statsB, 'data')
  statsB = statsB.data;
elseif options.keepDataWithNoStats
  return % statsB is missing, so keep origResults
else
  results = struct(); % statsB is missing, so call results unstable
  return
end

statFields = fieldnames(statsA);
assert( isequal(statFields, fieldnames(statsB)), ...
        'Different stats files have different kinds of saved stats!' )
for n = 1:length(statFields)
  sName = statFields{n};
  if strcmpi(sName, 'Time')
    continue  % don't check if time has changed :)
  end
  maxChange = getMaxChange(sName, options);
  A = statsA.(sName); medianA = median(A); maxA = max(A) ; minA = min(A);
  B = statsB.(sName); medianB = median(B); maxB = max(B) ; minB = min(B);
  changeA = 2 * (maxA - minA) / (abs(maxA + minA) + 1e-15);
  changeB = 2 * (maxB - minB) / (abs(maxB + minB) + 1e-15);
  changeAB = 2 * abs(medianA - medianB) / (abs(medianA + medianB) + 1e-15);
  if changeA > maxChange || changeB > maxChange || changeAB > maxChange
    % too much change, so return an empty structure (nothing is reliable)
    results = struct();
    return
  end
end
% these results passed the test, so return them as reliable
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxChange = getMaxChange(statName, options)
statName = [upper(statName(1)), lower(statName(2:end))];
maxChange = options.(sprintf('maxChange_%s', statName));
end
