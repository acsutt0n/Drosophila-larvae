% PlotStatsGraphs(ephysData, calciumCondition, holdingCondition, options)
% Make a plot from the statistics files.
%  INPUTS:
%  -ephysData: list of structs with FI information. Output of GetEphysData
%  -calciumCondition: string, 'Calcium' or 'No_Calcium'
%  -stat: string, 'Cm', 'Rm', 'Ra'
%  OPTIONAL INPUTS:
%  -colormap: string, name of matlab colormap (default='lines')
%  -title: string, title of plot (default='') if empty, a sensible plot
%          title will be generated from calcium and holding conditions
%  -marker: string, marker for FI data (default='-')
%  OUTPUTS:
%  -statistic: numStatsxnumData array of the values of the chosen
%              statistic
%  -independentValues: 1xnumData cell array of  strings,independent
%                      variable values (e.g. genotype of each cell)
function varargout = PlotStatsGraphs_CsChrimson( ephysData,...
                                         calciumCondition, stat, varargin )
  % define passable options and default values
  parser = inputParser();
  parser.addParameter('title', '')
  % set this to 'preparation_condition'
  parser.addParameter('independentVar', 'genotype')
  parser.addParameter('HoldingCondition', 'No_Holding')
  parser.addParameter('plotFunction', 'scatterPlot')
  
  % get the options
  parser.parse(varargin{:})
  options = parser.Results;
  
  iVar = options.independentVar;
  
  if isempty(options.title)
    options.title = [stat, ' - ', calciumCondition];
  end
  
  % create figure
  fig = NamedFigure(options.title) ; fig.WindowStyle = 'docked';
  clf( fig )
  
  independentValues = { ephysData.(iVar) };
  holdingCondition = options.HoldingCondition;
  
  % Find stats that aren't recorded by the atf files:
  if strcmp(stat, 'vHolding' )
      location = { calciumCondition, holdingCondition, 'ISteps', 'vHolding' };
      units = 'mV';
  elseif strcmp( stat, 'rheobase' )
      location = { calciumCondition, holdingCondition, 'Ramp', 'rheobase' };
      unitsLoc = { calciumCondition, holdingCondition, 'Ramp', 'rheobaseUnits' };
  elseif strcmpi( calciumCondition, 'Calcium' )
      location = { calciumCondition, 'stat_file_1', 'data', stat, @median };
      unitsLoc = { calciumCondition, 'stat_file_1', 'units', stat };
      afterLocation = { calciumCondition, 'stat_file_2', 'data', stat, @median };
  else
      location = { calciumCondition, 'stat_file_1', 'data', stat, @median };
      unitsLoc = { calciumCondition, 'stat_file_1', 'units', stat };
      afterLocation = { calciumCondition, 'stat_file_2', 'data', stat, @median };
  end
  
  % get units from unitsLoc if needed
  if ~exist( 'units', 'var' ) || isempty( units )
    units = arrayfun( @(d) GetDataFromLocation( d, unitsLoc ), ...
                      ephysData, 'UniformOutput', false );
    ind = find( cellfun( @(u) ischar(u), units ), 1 );
    units = units{ind};
  end

  % fix underscores in title
  options.title = strrep( options.title, '_', '\_' );

  %{
  % primer on how arrayfun works below:
  myFunc = @(d) GetDataFromLocation( d, location );
  % this code doesn't need to know how myFunc works:
  statistic = zeros( size( ephysData ) );
  for n = 1:numel( ephysData )
    d = ephysData(n);
    %statistic(n) = GetDataFromLocation( d, location );
    statistic(n) = myFunc( d );
  end
  %}
  
  
  
  if exist( 'afterLocation', 'var' )
    % statistic taken twice, plot before after, and change
    beforeStatistic = arrayfun( @(d) GetDataFromLocation( d, location ), ...
                                ephysData );
    afterStatistic = arrayfun( @(d) GetDataFromLocation( d, afterLocation ), ...
                               ephysData );
    statistic = [ beforeStatistic ; afterStatistic ];
    deltaStatistic = (afterStatistic - beforeStatistic) ...
                     .* 2 ./ abs( afterStatistic + beforeStatistic );
    % call specified plot function to plot before and after
    ax1 = subplot( 3,1,1:2, 'Parent', fig ) ; hold( ax1, 'on' )
    feval( options.plotFunction, ax1, beforeStatistic, independentValues, ...
           'color', 'b' );
    feval( options.plotFunction, ax1, afterStatistic, independentValues, ...
           'color', 'r' );
    legend( ax1, { 'before', 'after' }, 'Location', 'Best' )
    
    % plot change in statistic
    ax2 = subplot( 3,1,3, 'Parent', fig );
    feval( options.plotFunction, ax2, deltaStatistic, independentValues, ...
           'color', 'r' );
    ylabel( ax2, 'Fold change' )
    
  else % statistic taken just once, plot it
    statistic = arrayfun( @(d) GetDataFromLocation( d, location ), ...
                          ephysData );
  
    % call specified plot function
    ax1 = subplot( 1,1,1, 'Parent', fig );
    feval( options.plotFunction, ax1, statistic, independentValues );
  end
  
  ylabel( ax1, units )
  title( ax1, options.title )
  
  if nargout == 0
    varargout = {};
  else
    varargout = {statistic, independentValues};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make scatter plot for y-data with different independent variable values
function scatterPlot( varargin ) %#ok<DEFNU>
  if isa( varargin{1}, 'matlab.graphics.axis.Axes' )
    [ax, y, iVar] = varargin{1:3};
    varargin(1:3) = [];
  else
    ax = gca();
    [y, iVar] = varargin{1:2};
    varargin(1:2) = [];
  end
  parser = inputParser();
  parser.addParameter('marker', '')
  parser.addParameter('color', 'b')
  parser.parse( varargin{:} )
  options = parser.Results;
  
  if isempty( options.marker )
    markerOrder = { 'o', '^', 's', 'v', 'p', 'h', '<', '>' };
    markerInd = mod( numel( ax.Children ), numel( markerOrder ) ) + 1;
    options.marker = markerOrder{markerInd};
  end
  
  % get all the independent variable values and corresponding indices
  [allIVar, ~, iInd] = unique( iVar );
  numIVar = numel( allIVar );
  xRnd = 0.05 .* 2.0 .* ( rand( size(iInd) ) - 0.5 );
  plot( ax, iInd + xRnd, y, options.marker, 'Color', options.color )
  ax.XTick = 1:numIVar;
  ax.XTickLabel = allIVar;
  ax.TickLabelInterpreter = 'none';
  xlim( ax, [0.5 numIVar + 0.5] )
end