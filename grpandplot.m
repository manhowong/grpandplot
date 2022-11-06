function tile = grpandplot(data,yCol,varargin)
%% Group data automatically and plot scatter, box plot and/or violin plot.
% This function allows data grouping by multiple factors, supporting up to
% three-way grouping. Data can be grouped into groups along x-axis, into 
% different colors and/or into separate tiles (i.e. axes) of the same 
% figure. Different combination of data objects (i.e. scatter (data points)
% , box plot, violin plot, n number) can be plotted on top of each other.
% Figures are configurable with MATLAB's graphics methods and properties.

% Man Ho Wong (2022).
% -------------------------------------------------------------------------
% Input: - data [Table]
%               contains both data to be plotted (Y values) and
%               grouping factors (if any). Y values and grouping factors
%               are organized in separate columns, and they must have the
%               same length. If no grouping factor is provided, data will
%               be plotted as one group.
%        - yCol [String (scalar) | Character array | Positive integer]
%               specifies the name (or the index number) of the column
%               containing the Y values in the table 'data'.
%        ------------------------------------------------------------------
%        Optional name-value pair input:
%
%        - parent [graphics object]
%                 The graphics object (an axes object, a figure or a tiled
%                 chart layout) where the tiles(s) will be drawn in. This
%                 object must contain sufficient space. For example, if
%                 2-by-2 tiles (i.e. four axes objects) will be drawn, a 
%                 2-by-2 tile space must be available in 'parent'.
%                 - If 'parent' is not specified, a new figure with be
%                   created. The tile layout (number of rows and columns)
%                   will follow 'nRow' and 'nCol' if they are provided.
%                 - If 'parent' is a figure object, the tile layout will
%                   follow 'nRow' and 'nCol' if they are provided.
%                 - If 'parent' is a tiled chart layout, the tile layout
%                   will follow the layout of 'parent', ignoring 'nRow' and
%                   'nCol' settings.                 
%                 - If 'parent' is an axes object, grouping data in tiles
%                   will not be supported.
%        - nRow [Positive integer]
%               Number of tile rows (if grouping data into tiles)
%        - nCol [Positive integer]
%               Number of tile columns (if grouping data into tiles)
%        - xFactor [String (scalar) | Character array] 
%                  Factor for grouping data into groups along x-axis
%        - cFactor [String (scalar) | Character array] 
%                  Factor for grouping data into separate colors
%        - tFactor [String (scalar) | Character array] 
%                  Factor for grouping data into separate tiles
%        - xOrder  [String | Character array | Cells of strings]
%                  Custom order of 'xFactor' levels
%        - cOrder  [String | Character array | Cells of strings]
%                  Custom order of 'cFactor' levels
%        - tOrder  [String | Character array | Cells of strings]
%                  Custom order of 'tFactor' levels
%        - yTitle [String (scalar) | Character array]
%                 By default, y-axis title is the 'yCol' input. Use this
%                 option to overwrite the default y-axis title.
%        - xAxisMode [0, 1 or 2] (default: 0)
%                    Set 'xAxisMode' to 1 to turn on multi-level x-axis
%                    labeling. X-axis will be labeled by color grouping
%                    factors first and then by x-axis grouping factors.
%                    Option 2 is similar to 1, except that each label is
%                    an independent object which allows customization.
% 
%        The following options control what data objects to include.
%        [logical (accept true, false, 0, 1)]
%        - showPnt : Data points (default: true)
%        - showBox : Box plot (default: true)
%        - showVln : Violin plot (default: false)
%        - showOutlier : Box plot outliers (default: false)
%        - showXLine : Vertical line separating groups allow x-axis. 
%                      Only works if grouping by xFactor is active.
%                      (default: false)
%        - showLegend : Legend (default: true)
%        - showNum : n number (number of data points) (default: false)
%
%        Position and size:
%        - numYPos [numeric]
%                  Y-position of n number.
%        - pntOnTop [logical (accept true, false, 0, 1)] (default: true)
%                   Show data points on top of other graphics objects.
%        - pntSize [positive number] (default: 20)
%                  Size of data point marker. 
%        - w [positive number] (default: 0.7 for 0 or 1 color group)
%            width of data object; relative to width between 2 x-axis labels
%        - gap [positive number] (default: 2.2 for 0 or 1 color group)
%              space between data objects of same x-axis groups if 'cFactor'
%              is given; relative to number of color groups
%        - xSpace [positive number] (default: 0)
%                 space between data objects and ends of x-axis; use this
%                 to adjust x-axis length (i.e. horizontal scale of figure)
% 
%        Color:
%        - cmap [Rows of RGB triplets] (default : MATLAB's "lines")
%               Color map with each color represented by an RGB triplet.
%               e.g. [0 0 0; 1 0 0; 0 0 1] = black, red, blue. 'cmap'
%               specifies the colors of the data objects if data is grouped
%               into different colors; Otherwise, the first color in 'cmap'
%               will be used for all data objects.
%        - pntFillC [RGB triplet]
%                   Overwrite fill color of data points with custom color.
%        - boxFillC [RGB triplet]
%                   Overwrite fill color of box plots with custom color.
%        - vlnFillC [RGB triplet]
%                   Overwrite fill color of violin plots with custom color.
%        - pntEdgeC [RGB triplet] (default: black)
%                   Overwrite edge color of data points with custom color.
%        - boxEdgeC [RGB triplet]
%                   Overwrite edge color of box plots with custom color.
%        - vlnEdgeC [RGB triplet]
%                   Overwrite edge color of violin plots with custom color.
%        - pntAlpha [Positive number between 0-1] (default: 1)
%                   Transparency of data points. 0 = full transparency.
%        - boxAlpha [Positive number between 0-1] (default: 0.5)
%                   Transparency of box plots.  0 = full transparency.
%        - vlnAlpha [Positive number between 0-1] (default: 0.5)
%                   Transparency of violin plots. 0 = full transparency.
% -------------------------------------------------------------------------
% Output: - new figure or Axes object(s) in the 'parent' object
%         - tile [Axes object or an array of Axes objects] (Optional)
%           Handle of an Axes object or an array of Axes objects in the 
%           'parent' object (if specified) or in a new figure. Use this
%           handle to configure an Axes object with MATLAB's graphics
%           methods and properties.
% Example: See grpandplot_doc.mlx.

%% 1. Input parsing

% Validation functions
validNum = @(x) isnumeric(x) && isscalar(x);
validPosNum = @(x) validNum(x) && (x>0);
validPosInt = @(x) validPosNum(x) && (mod(x,1)==0);
validPosFrc = @(x) validNum(x) && (x>=0) && (x<=1);
validLogical = @(x) islogical(x) || islogical(logical(x));
validName = @(x) isStringScalar(x) || ischar(x) || ~isempty(x) || validPosInt(x);
validStrs = @(x) isstring(x) || ischar(x) || iscellstr(x);
validColor = @(x) ~isempty(validatecolor(x));
validColors = @(x) ~isempty(validatecolor(x,'multiple'));

% Build input parser
p = inputParser;
addRequired(p,'data', @istable);
addRequired(p,'yCol', validName);
addParameter(p,'parent',[], @isobject);
addParameter(p,'nRow', 1, validPosInt);
addParameter(p,'nCol', 1, validPosInt);
addParameter(p,'xFactor', '', validName);
addParameter(p,'cFactor', '', validName);
addParameter(p,'tFactor', '', validName);
addParameter(p,'xOrder', [], validStrs);
addParameter(p,'cOrder', [], validStrs);
addParameter(p,'tOrder', [], validStrs);
addParameter(p,'yTitle', yCol, validName);
addParameter(p,'showPnt', true, validLogical);
addParameter(p,'showBox', true, validLogical);
addParameter(p,'showVln', false, validLogical);
addParameter(p,'showOutlier', false, validLogical);
addParameter(p,'showXLine', false, validLogical);
addParameter(p,'showLegend', true, validLogical);
addParameter(p,'showNum', false, validLogical);
addParameter(p,'numYPos', 0, validNum);
addParameter(p,'pntOnTop', true, validLogical);
addParameter(p,'pntSize', 20, validPosNum);
addParameter(p,'xAxisMode', 0, validNum);
addParameter(p,'w', [], validPosNum);
addParameter(p,'gap', [], validPosNum);
addParameter(p,'xSpace', 0, validPosNum);
addParameter(p,'cmap', colormap("lines"), validColors);
addParameter(p,'pntFillC', [], validColor);
addParameter(p,'boxFillC', [], validColor);
addParameter(p,'vlnFillC', [], validColor);
addParameter(p,'pntEdgeC', 'k', validColor);
addParameter(p,'boxEdgeC', [], validColor);
addParameter(p,'vlnEdgeC', [], validColor);
addParameter(p,'pntAlpha', 1, validPosFrc);
addParameter(p,'boxAlpha', 0.5, validPosFrc);
addParameter(p,'vlnAlpha', 0.5, validPosFrc);

% Parse input
parse(p,data,yCol,varargin{:});
data = p.Results.data;
yCol = p.Results.yCol;
parent = p.Results.parent;
nRow = p.Results.nRow;
nCol = p.Results.nCol;
yTitle = p.Results.yTitle;
xFactor = p.Results.xFactor;
cFactor = p.Results.cFactor;
tFactor = p.Results.tFactor;
xOrder = p.Results.xOrder;
cOrder = p.Results.cOrder;
tOrder = p.Results.tOrder;
showPnt = p.Results.showPnt;
showBox = p.Results.showBox;
showVln = p.Results.showVln;
showOutlier = p.Results.showOutlier;
showXLine = p.Results.showXLine;
showLegend = p.Results.showLegend;
showNum = p.Results.showNum;
numYPos = p.Results.numYPos;
pntOnTop = p.Results.pntOnTop;
pntSize = p.Results.pntSize;
xAxisMode = p.Results.xAxisMode;
w = p.Results.w;
gap = p.Results.gap;
xSpace = p.Results.xSpace;
cmap = p.Results.cmap;
pntFillC = p.Results.pntFillC;
boxFillC = p.Results.boxFillC;
vlnFillC = p.Results.vlnFillC;
pntEdgeC = p.Results.pntEdgeC;
boxEdgeC = p.Results.boxEdgeC;
vlnEdgeC = p.Results.vlnEdgeC;
pntAlpha = p.Results.pntAlpha;
boxAlpha = p.Results.boxAlpha;
vlnAlpha = p.Results.vlnAlpha;

%% 2. Data preprocessing

%--------------------------------------------------------------------------
% Convert factor levels to categories
factors = {xFactor, tFactor, cFactor};
orders = {xOrder, tOrder, cOrder};
for i = 1:3
    if ~isempty(factors{i})
        % Convert factor levels to cellstr then to categorical datatype
        data = convertvars(data,factors{i},'cellstr');
        data = convertvars(data,factors{i},'categorical');
        if ~isempty(orders{i})
            % Reorder factor levels as specified by user
            data.(factors{i}) = reordercats(data.(factors{i}),orders{i});
        end
    end
end

%--------------------------------------------------------------------------
% Get factor levels (i.e. categories) and number of levels
nXLvl = 1;  % If xFactor is not provided, put all datapoints in one group.
xLvls = 'X';  % If xFactor is not provided, put all datapoints in one group.
if ~isempty(xFactor)
    xLvls = categories(data.(xFactor));
    nXLvl = length(xLvls);
end
nTLvl = 1;  % default
if ~isempty(tFactor)
    tLvls = categories(data.(tFactor));
    nTLvl = length(tLvls);
end
nCLvl = 1;  % default
if ~isempty(cFactor)
    cLvls = categories(data.(cFactor));
    nCLvl = length(cLvls);
end

% if number of colors in cmap is less than nCLvl, repeat the colors in cmap
nColors = height(cmap);
if nColors < nCLvl
    m = ceil(nCLvl/nColors);
    cmap = repmat(cmap, m,1);
    warning(['Number of colors in cmap (%d) is less than number of color ' ...
        'factor levels (%d). Colors are repeated.'], length(cmap), nCLvl);
end

%--------------------------------------------------------------------------
% Group data according to grouping factors
if nTLvl > 1
    for i = 1:nTLvl
        % Data grouped by tFactor
        dataByTile{i} = data(data.(tFactor)==tLvls(i),:);
        if nCLvl > 1
            for j = 1:nCLvl
                % Data grouped by tFactor and then by cFactor
                dataByTileByColor{i,j} = ...
                    dataByTile{i}(dataByTile{i}.(cFactor)==cLvls(j),:);
                    % Note: Plotting functions will group data by xFactor
            end
        end
    end
elseif nCLvl > 1
    for j = 1:nCLvl
        % Data grouped by cFactor
        dataByColor{1,j} = data(data.(cFactor)==cLvls(j),:);
    end
end

%--------------------------------------------------------------------------
% Assign grouped data to ydata for plotting
if nTLvl > 1
    ydata = dataByTile';  % transpose data to match indexing style later
    if nCLvl > 1
        ydata = dataByTileByColor;
    end
elseif nCLvl > 1
    ydata = dataByColor;
else
    ydata{1} = data; % assign data as a cell to match indexing style later 
end

%% 3. Plotting

% Get number of tiles
nTiles = nRow*nCol;

% If there are more tile levels than number of tiles set by user (through
% nRow and nCol), reset nRow and nCol so that there will be enough tiles
if nTLvl > nTiles
    nRow = 1;
    nCol = nTLvl;
    nTiles = nRow*nCol;
    warning(['nRow and nCol not provided, or number of tile levels exceeds number of ' ...
        'tiles (nRow*nCol): Tile layout are automatically set to fit all tiles.']);
end

% if no parent is specified, a new figure will be created.
if isempty(parent)
    parent=figure;
end

% if there's only layout but no figure, a new figure will be created.
if isgraphics(parent,'tiledlayout') && parent.Parent.Tag == "EmbeddedFigure_Internal"
    greatparent=figure;
    parent.Parent = greatparent;
end

if isgraphics(parent,'figure')
    parent = tiledlayout(nRow,nCol,'Padding','loose','TileSpacing','loose');
    for i = 1:nTiles
        tile(i) = nexttile(parent);
    end
elseif isgraphics(parent,'tiledlayout')
    for i = 1:nTiles
        tile(i) = nexttile(parent);
    end
elseif isgraphics(parent,'axes') && nTLvl == 1
    tile(1) = parent;
else
    warning(['The provided ''Parent'' object does not support tiled layout. ' ...
             'Please provide another object as the ''Parent'', ' ...
             'or disable grouping by tile factors.']);
    return;
end

% Set data object width and spacing
if isempty(w)
    w = 0.7/nCLvl;  % adjust width automatically by nCLvl
end
if isempty(gap)
    gap = 2.2/nCLvl;  % adjust gap automatically by nCLvl
end

% Plot each tile one by one -----------------------------------------------
for i = 1:nTLvl

    % Set current Axes to tile(i) (so it is the target Axes for plotting)
    %  *One can also specify the target Axes in individual plot functions
    %  (eg boxchart, swarnchart and text), except violin which does not
    %  have this option.
    axes(tile(i));
    hold(tile(i),"on")

    % x-position offset for each color group (so that they don't overlap)
    xOffset = -(gap*0.2)*(nCLvl-1) : 2*(gap*0.2) : (gap*0.2)*(nCLvl-1);

    % For each tile, plot each color group one by one
    for j = 1:nCLvl

        % Get x-position of each group of data on the plot by converting
        %   x factor levels to numbers and use them as the x-positions
        xIDs = [];
        if ~isempty(xFactor)
            for xx = 1:nXLvl
                idx = find(ydata{i,j}.(xFactor)==xLvls(xx));
                xIDs(idx,1) = xx;
            end
        else
            % if no grouping by x, assign all data points to same group no.
            xIDs = ones(height(ydata{i,j}),1);
        end

        xPos = xIDs + xOffset(j);  % a list of x-positions
        y = ydata{i,j}.(yCol);  % data for plotting data points and box plot

        % Show raw data points at bottom layer-----------------------------
        if showPnt == true && pntOnTop == false
            pntObj(i,j) = swarmchart(xPos, y, MarkerFaceColor=cmap(j,:), ...
                                     MarkerEdgeColor=pntEdgeC, ...
                                     MarkerFaceAlpha=pntAlpha);
            % Configure pntObj
            pntObj(i,j).XJitterWidth= w;
            pntObj(i,j).SizeData = pntSize;
            if ~isempty(pntFillC), pntObj(i,j).MarkerFaceColor = pntFillC; end
        end

        % Show box plot----------------------------------------------------
        if showBox == true
            boxObj(i,j) = boxchart(xPos, y, BoxFaceColor=cmap(j,:), ...
                                   BoxFaceColorMode="manual");
            % Configure plotObj
            if showOutlier == false                
                boxObj(i,j).MarkerStyle = 'none';  % hide outlier points
            end
            boxObj(i,j).MarkerColor = 'r';
            boxObj(i,j).BoxWidth = w;
            boxObj(i,j).BoxFaceAlpha = boxAlpha;
            boxObj(i,j).BoxLineColor = cmap(j,:);
            boxObj(i,j).WhiskerLineColor = cmap(j,:);
            if ~isempty(boxFillC), boxObj(i,j).BoxFaceColor = boxFillC; end
            if ~isempty(boxEdgeC)
                boxObj(i,j).BoxLineColor = boxEdgeC;
                boxObj(i,j).WhiskerLineColor = boxEdgeC;
            end
        end

        % Show violin plot-------------------------------------------------
        if showVln == true
            for xx = 1:nXLvl
                % Plot each x level separately
                % (Unlike boxchart/swarmchart, automatic grouping and
                %  plotting by x level is not available for violin)
                if ~isempty(xFactor)
                    xLvlFilter = ydata{i,j}.(xFactor)==xLvls(xx);
                    yViolin = ydata{i,j}.(yCol)(xLvlFilter);
                else
                    yViolin = ydata{i,j}.(yCol)(:);
                end                
                if ~isempty(yViolin)  % don't plot if no data for current x level
                    vlnObj(i,j) = violin(xx + xOffset(j), yViolin, linecolor=cmap(j,:),...
                                         scaling=w, withmdn=1, cutoff=0);
                    vlnObj(i,j).FaceColor = cmap(j,:);
                    vlnObj(i,j).EdgeColor = cmap(j,:);
                    vlnObj(i,j).FaceAlpha = vlnAlpha;
                    if ~isempty(vlnFillC), vlnObj(i,j).FaceColor = vlnFillC; end
                    if ~isempty(vlnEdgeC), vlnObj(i,j).EdgeColor = vlnEdgeC; end
                end
            end
        end

        % Show raw data points on top of other objects---------------------
        if showPnt == true && pntOnTop == true
            pntObj(i,j) = swarmchart(xPos, y, MarkerFaceColor=cmap(j,:), ...
                                     MarkerEdgeColor=pntEdgeC, ...
                                     MarkerFaceAlpha=pntAlpha);
            % Configure pntObj
            pntObj(i,j).XJitterWidth= w;
            pntObj(i,j).SizeData = pntSize;
            if ~isempty(pntFillC), pntObj(i,j).MarkerFaceColor = pntFillC; end
        end

        % Show n number----------------------------------------------------
        if showNum == true
            for xx = 1:nXLvl
                % Add n number for each x level separately
                % (Unlike boxchart/swarmchart, automatic grouping and
                % plotting by x level is not available for text annotation)
                xLvlFilter = ydata{i,j}.(xFactor)==xLvls(xx);
                yAtX = ydata{i,j}.(yCol)(xLvlFilter);
                n = num2str( length( yAtX(~isnan(yAtX)) ) );
                % n = num2str(length(xIDs(xIDs==xx)));
                if n ~= 0   % don't show n number if n = 0
                    text(xx + xOffset(j), numYPos, n, ...
                        FontSize=8.55, HorizontalAlignment="center");
                end
            end
        end
    end

    if showXLine == true
        xline(tile(i),1.5:1:nXLvl,'LineStyle',':');
    end
end

%% 4. Figure settings
% Specify the tile number to configure in the functions below; to configure
% all tiles, enter tile(:) (only allowed in some functions) or loop tiles
% one by one.

% Tile titles -------------------------------------------------------------
for i = 1:nTLvl
    if exist('tLvls','var')
        title(tile(i), tFactor + ": " + tLvls(i));
    end
end

% Y-axis-------------------------------------------------------------------
ylabel(tile(:), yTitle);  % y-axis title
% tile(end).YLabel.Visible = false;  % remove y-axis title
ylim(tile(:), [0 max(data.(yCol))*1.1]);  % y-axis starts at 0, ends 10% more than max. value
linkaxes(tile(:), 'y');  % link y-axes of all tiles in tile set


% X-axis scale ------------------------------------------------------------
% adjust x-axis length via x-axis range
if nCLvl > 1
    xlim(tile(:), [1 + min(xOffset) - 1.2*w - xSpace, ...
               nXLvl + max(xOffset) + 1.2*w + xSpace]);
else
    xlim(tile(:), [1 - w - xSpace, nXLvl + w + xSpace]);
end

% Legend-------------------------------------------------------------------
% Decide which data objects are used in legend
%  Priority: scatter > box plot > violin plot

if nCLvl > 1 
    if exist("pntObj", "var") && isempty(pntFillC)
        obj = pntObj(1,:);
    elseif exist("boxObj", "var") && (isempty(boxEdgeC) || isempty(vlnFillC))
        obj = boxObj(1,:);
    elseif exist("vlnObj","var") && (isempty(vlnEdgeC) || isempty(vlnFillC))
        obj = vlnObj(1,:);
    end
end

% Add legend for obj to last tile; use cLvls as the text labels
if exist("obj","var") && showLegend == true
    legend(tile(end),obj,cLvls,Location='bestoutside',Box='off');
end

% Other adjustments -------------------------------------------------------

for t = 1:nTiles
    set(tile(t).XAxis,'TickDir','out');  % x-tick direction
    set(tile(t).YAxis,'TickDir','out');  % y-tick direction
    set(tile(t),'TickLength', [0.02, 0.02]);
    tile(t).FontSize = 11;
    tile(t).FontName = 'Arial';
    tile(t).LineWidth = 1;
end

% X-axis labels -----------------------------------------------------------
% Configure x-axis labeling last after everything else is set since
% positioning of labels here is normalized to the figure's dimensions


if xAxisMode == 0 || nCLvl == 1 % default mode
    % Show ticks only at xLvls
    xticks(tile(:), 1:1:nXLvl);   % x-axis tick positions
    xticklabels(tile(:), xLvls);  % use xLvls as tick labels 

elseif xAxisMode == 1 && nCLvl > 1  % two-level x-axis labeling
    % get x-axis tick positions from xOffset
    allXPos = [];
    for xx = 1:nXLvl
        xPos = xx + xOffset;
        allXPos = [allXPos xPos];
    end
    xticks(tile(:), allXPos);

    % 1st-level labels (groups by cFactor)
    labels = repmat(cLvls,[1,nXLvl]);
    xticklabels(tile(:), labels);

    % add 2nd-level labels and group lines
    for t = 1:length(tile)
        h = tile(t);
        % add an empty x-axis title to make space and use its y-position as
        % the position for 2nd-level label group lines
        fakeTitles(t) = xlabel(h,{' ' ' '});
        grpLineYPos = fakeTitles(1).Position(2)*1.2;
        xLabelYPos = grpLineYPos*1.2;
        % add 2nd-level label group lines
        for i = 1:nCLvl:length(allXPos)
            line(h,[allXPos(i)-w/2 allXPos(i+nCLvl-1)+w/2], ...
                 [grpLineYPos grpLineYPos], ...
                 LineWidth=1, Color='k', Clipping='off');
        end
        % add 2nd-level labels (groups by xFactor)
        for xx = 1:nXLvl
            text(h, xx, xLabelYPos, xLvls(xx), FontSize=11, ...
                 HorizontalAlignment="center", Clipping='off');
        end
    end

elseif xAxisMode == 2 && nCLvl > 1  % two-level x-axis labeling (customizable)
    % get x-axis tick positions from xOffset
    allXPos = [];
    for xx = 1:nXLvl
        xPos = xx + xOffset;
        allXPos = [allXPos xPos];
    end
    xticks(tile(:), allXPos);
    % remove default x tick labels
    xticklabels(tile(:), '');
    % add labels as text objects (each label is customizable)
    for t = 1:length(tile)
        h = tile(t);
        yHeight = ylim;
        % 1st-level labels (groups by cFactor)
        labels = repmat(cLvls,[1,nXLvl]);
        for i = 1:length(allXPos)
            text(h, allXPos(i), -yHeight(2)/8, labels(i), ...
                 FontSize=9.1, HorizontalAlignment="center");
        end 
        % add 2nd-level labels (groups by xFactor)
        for xx = 1:nXLvl
            text(h, xx, -yHeight(2)*0.34, upper(xLvls(xx)), FontSize=11, ...
                 HorizontalAlignment="center", Clipping='off');
        end
        % add 2nd-level label group lines
        for i = 1:nCLvl:length(allXPos)
            line(h,[allXPos(i)-w/2 allXPos(i+nCLvl-1)+w/2], ...
                [-yHeight(2)*0.23 -yHeight(2)*0.23], ...
                Clipping='off',LineWidth=1,Color='k');
        end
    end
end