function [outDet, labels] = brushDet(DATA, inDet, varargin)
% [outDet, labels] = brushDet(DATA, inDet, paramFile) runs GUI which 
% allows user to delete and label detections on top of time series data
% -DATA{i} is the time series data for the ith receiver
% -DATA{i}(1,:) is the time vector for the ith receiver
% -DATA{i}(2,:) is the amplitude vector for the ith receiver
% -inDet{i} are the initial detections for the ith receiver
% -inDet{i}(1,:) are the detection times
% -(optional) paramFile is the path to a parameter file which allows various
% settings for the interface
% -inDet{i}(2,:) are the detection amplitudes
% -outDet are the retained detections, in the same format as inDet
% -labels{i} are the labels assigned to detections

global brushing

brushing.DATA = DATA;           % acoustic data
brushing.nAxes = numel(DATA);   % number of axes/receivers to display
brushing.inDet = inDet;         % initial detections
brushing.outDet = inDet;        % retained detections after adjustments/brushing

% determine if param file was specified:
if nargin < 2
    errordlg('Not enough input arguments', 'Error')
elseif nargin == 2 % no parameter file specified
    % check for existing default param file
    if isfile('brushDet.params')
        loadParams('brushDet.params')
    else
        makeParamFile
    end
else nargin == 3 % parameter file specified
    
end

% test whether label color legend exists
figCol = findall(0, 'Type', 'figure', 'name', 'Legend of Label Colors');
if isempty(figCol) 
    generateColorSchemeLegend(brushing) % if legend doesn't exist, generate it
end

% test whether Keystroke Command Legend exists
figKey = findall(0, 'Type', 'figure', 'name', 'Keystroke Command Legend');
if isempty(figKey)
    generateKeystrokeLegend(brushing) % if legend doesn't exist, generate it
end

run_brushDet % fun the GUI

qselect = input('\nEnter ''q'' to quit: ', 's');
outDet = brushing.outDet;
labels = brushing.labels;

end


function run_brushDet
% runs the acutal brushDet GUI

% initialize global values
global brushing

brushing.tRange = [nan, nan];        % limits of x axis on plot
for i = 1:brushing.nAxes
    brushing.colorIndex{i} = 2.*ones(length(brushing.inDet{i}), 1); % index of brushing.params.colorMat row for unlabeled detections
    brushing.tooLow{i} = [];    % points removed due to being lower than slider threshold value

    brushing.tRange(1) = min([min(brushing.DATA{i}(1,:)), brushing.tRange(1)]); % minimum time on plot
    brushing.tRange(2) = max([max(brushing.DATA{i}(1,:)), brushing.tRange(2)]); % maximum time on plot

    brushing.labels{i} = num2str(zeros(length(brushing.inDet{i}), 1)); % initialize labels
end


uif = findall(0, 'Type', 'figure', 'name', 'Brush Detections'); % activate Brush Detections GUI (if it exists)

if isempty(uif) % If Brush Detections doesn't exist, run init_brushDet
    [uif, ax, sld] = init_brushDet;

else
    if brushing.nAxes~=numel(uif.Children.Children) % detect a change in the number of receivers/axes
        close(uif)
        [uif, ax, sld] = init_brushDet;
    end
    for i = 1:brushing.nAxes
        ax(i) = uif.Children.Children(i).Children(2);
        sld(i) = uif.Children.Children(i).Children(1);
    end
end

for i = 1:brushing.nAxes
    % set slider callback function
    set(sld(i), 'ValueChangedFcn', {@changeThresh, sld(i)})

    plot(ax(i), brushing.DATA{i}(1,:), brushing.DATA{i}(2,:), 'HandleVisibility', 'off', 'color', brushing.params.colorMat(1,:))

    if ~isempty(brushing.inDet{i}) % determine if there are detections on this receiver
        hold(ax(i), 'on')
        scatter(ax(i), brushing.outDet{i}(1, :), brushing.outDet{i}(2, :), 100,...
            brushing.params.colorMat(brushing.colorIndex{i}, :), 'x', 'linewidth', 2) % plot detections
        hold(ax(i), 'off')
        xlim(ax, brushing.tRange)

        if length(brushing.outDet{i})>1
            % adjust slider range according to detection amplitudes
            uif.Children.Children(i).Children(1).Limits = [min(brushing.outDet{i}(2,:)), max(brushing.outDet{i}(2,:))];

            % Set slider value to minimum detection amplitude
            uif.Children.Children(i).Children(1).Value = min(brushing.outDet{i}(2,:));

            b = brush(uif);
            b.Color = 'g';
            b.Enable = 'on';
        else
            fprintf('\ninsufficient detections on Receiver %d\n', i)
        end
    else
        fprintf('\nno detections on Receiver %d\n', i)
    end
end

hManager = uigetmodemanager(uif);
[hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2 (on 2014b or later)
set(uif, 'KeyPressFcn', @keyPressCallback);

end

function [uif, ax, sld] = init_brushDet
% Initializes the Brush Detections plot.
% nAxes = number of axes in the plot (number of receivers)
% uif is the figure handle
% ax(i) is the axes handle for receiver i
global brushing
uif = uifigure('Name', 'Brush Detections');
uif.Position = brushing.params.figPosition;
set(uif, 'ToolBar', 'none')
set(uif, 'MenuBar', 'none')
grd = uigridlayout(uif, [brushing.nAxes, 1]); % divide uifigure into nAxes X 1 grids

for i = 1:brushing.nAxes

    uipan(i) = uipanel(grd, 'Title', ['Receiver ', num2str(i)]); % create panel on each grid

end
drawnow
pause(1)
for i = 1:brushing.nAxes
    ax(i) = axes(uipan(i), 'Position', [0.08, 0.11, 0.88, 0.815]); % make plot axes on each panel

    % Set slider position:
    sldPos(1) = 10; % Distance from left edge of panel (pixels)
    sldPos(2) = 10;  % Distance from bottom edge of panel (pixels)
    sldPos(3) = 3;
    sldPos(4) = uipan(i).InnerPosition(4)-20;

    sld(i) = uislider(uipan(i), 'Orientation', 'vertical', ...
        'Position', sldPos);

end

end

function changeThresh(hObject, sldVal, handles)
% Detects a slider change and updates the detections
% sldVal.Value is the current slider value
% sldVal.PreviousValue is the previous slider value

global brushing

ax = hObject.Parent; % get the axis handle
iRec = str2double(ax.Title(end)); % Receiver number

if sldVal.Value>sldVal.PreviousValue % slider value increased
    
    Irem = find(brushing.outDet{iRec}(2,:)<sldVal.Value); % indices of detections which are lower than threshold
    
    brushing.tooLow{iRec} = [brushing.tooLow{iRec}, brushing.outDet{iRec}(:, Irem)]; % add points to be removed to brushing.tooLow
    % ^ I do this so the data are still saved in case the user lowers the
    % threshold again later
    
    brushing.outDet{iRec}(:, Irem) = [];
    brushing.labels{iRec}(Irem) = [];
    brushing.colorIndex{iRec}(Irem) = []; % remove color index of selected data
    
    set(ax.Children(2).Children, 'xdata', brushing.outDet{iRec}(1,:), 'ydata', brushing.outDet{iRec}(2,:), ...
        'cdata', brushing.params.colorMat(brushing.colorIndex{iRec}, :));

elseif sldVal.Value<sldVal.PreviousValue % slider value decreased
    Iadd = find(brushing.tooLow{iRec}(2,:)>=sldVal.Value); % indices of points to be added back onto plot
    
    % add those points back into plotted data
    brushing.outDet{iRec} = [brushing.outDet{iRec}, brushing.tooLow{iRec}(:,Iadd)]; % add removed points back into data
    brushing.labels{iRec} = [brushing.labels{iRec}; num2str(zeros(length(Iadd),1))];   % points being added in are unlabeled
    brushing.colorIndex{iRec} = [brushing.colorIndex{iRec}; 2.*ones(length(Iadd), 1)];   % remove color index of selected data

    % remove these points from brushing.tooLow
    brushing.tooLow{iRec}(:,Iadd) = [];

    set(ax.Children(2).Children, 'xdata', brushing.outDet{iRec}(1,:), 'ydata', brushing.outDet{iRec}(2,:), ...
        'cdata', brushing.params.colorMat(brushing.colorIndex{iRec}, :));

else % slider value did not change-- if you're here something went wrong
    errBox = msgbox('error: slider value didn''t change. I don''t know how you got here.', 'Magic Error');
end



end

%% kepPressCallback Function
function keyPressCallback(source,eventdata)
% runs when use presses a 
global brushing

% Receive keyboard input
key = eventdata.Key; % pressed key
% fprintf(['\n''', key, ''' key pressed\n']) % print message (for debugging)

for i = 1:brushing.nAxes % get highlighted points on each receiver
    selectedData = source.Children.Children(i).Children(2).Children.BrushData;
    Ind{i} = find(selectedData~=0); % indices of selected data
end

numkey = str2double(key); % convert keyboard input to a number (returns NaN if value is not a number)

if ~(isempty(numkey)||isnan(numkey))
    if numkey>8
        errBox = msgbox('error: whale number too high\nSelect a number between 1 and 8', 'Error');
    else
        for i = 1:brushing.nAxes
            if ~isempty(Ind{i})
                source.Children.Children(i).Children(2).Children.BrushData(Ind{i}) = []; % unselect all brushed points
                brushing.labels{i}(Ind{i}) = key; % set labels as the whalenumber
                brushing.colorIndex{i}(Ind{i}) = numkey+2; % set color of datapoints

                set(source.Children.Children(i).Children(2).Children, ...
                    'cdata', brushing.params.colorMat(brushing.colorIndex{i}, :));
            end
        end
    end
else

    switch key
        case 'd' % delete selected points
            for i = 1:numel(Ind)
                if ~isempty(Ind{i})

                    source.Children.Children(i).Children(2).Children.BrushData(Ind{i}) = []; % unselect all brushed points
                    brushing.outDet{i}(:,Ind{i}) = []; % remove selected data from dataset
                    brushing.labels{i}(Ind{i}) = [];   % remove labels of selected data
                    brushing.colorIndex{i}(Ind{i}) = []; % remove color index of selected data

                    set(source.Children.Children(i).Children(2).Children, ...
                        'xdata', brushing.outDet{i}(1,:), 'ydata', brushing.outDet{i}(2,:), ...
                        'cdata', brushing.params.colorMat(brushing.colorIndex{i}, :));

                end
            end
        case 'b' % label selected points as a buzz
            for i = 1:numel(Ind)
                if ~isempty(Ind{i})

                    source.Children.Children(i).Children(2).Children.BrushData(Ind{i}) = []; % unselect all brushed points
                    brushing.labels{i}(Ind{i}) = 'b'; % set labels as 'b'
                    brushing.colorIndex{i}(Ind{i}) = 11; % set color of datapoints
                    
                    set(source.Children.Children(i).Children(2).Children, ...
                        'cdata', brushing.params.colorMat(brushing.colorIndex{i}, :));

                end
            end
        case 'z' % turn on zoom
            
            % enable zoom functionality
            z = zoom(source); % get zoom object of figure

            z.Enable = 'on'; % turn on zoom
            z.Direction = 'in'; % set zoom to "in"
            
            % reenable keyPress (zoom and brush automatically disable this)
            hManager = uigetmodemanager(source);
            [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2 (on 2014b or later)
            set(source, 'KeyPressFcn', @keyPressCallback);

            
        case 'x' % turn off zoom
            z = zoom(source); % get zoom object of figure
            z.Enable = 'off';  % turn on zoom
            
            % reenable brush:
            b = brush(source);
            b.Color = 'g';
            b.Enable = 'on';

            % reenable keyPress (zoom and brush automatically disable this)
            hManager = uigetmodemanager(source);
            [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2 (on 2014b or later)
            set(source, 'KeyPressFcn', @keyPressCallback);
        case 'n' % select a new detection
            
            iRecCell = inputdlg('Which receiver?', 'Select New Detections');
            iRec = str2double(iRecCell{1});
            if iRec <= numel(source.Children.Children) && iRec~=0


                detPoint = drawpoint(source.Children.Children(iRec).Children(2)); % point struct of selected datapoint(s)
                
                
                distToDatapoints = sum((brushing.DATA{iRec} - detPoint.Position.').^2); % distance between selected point and all points in data
                [~, Imin] = min(distToDatapoints);

                newDet = brushing.DATA{iRec}(:, Imin);

                % add point into plotted data
                brushing.outDet{iRec} = [brushing.outDet{iRec}, newDet]; % add removed points back into data
                brushing.labels{iRec} = [brushing.labels{iRec}; num2str(0)];   % points being added in are unlabeled
                brushing.colorIndex{iRec} = [brushing.colorIndex{iRec}; 2];   % remove color index of selected data
                
                set(ax.Children(2).Children, 'xdata', brushing.outDet{iRec}(1,:), 'ydata', brushing.outDet{iRec}(2,:), ...
                    'cdata', brushing.params.colorMat(brushing.colorIndex{iRec}, :));

            else
                errBox = msgbox(['error: invalid receiver number. \nEnter a number from 1 to ', num2str(numel(source.Children.Children))], 'Error');
            end

            
            
    end
end
end

%%
function makeParamFile
% make a parameter file for function settings

global brushing

fid = fopen('brushDet.params', 'w');
fprintf(fid, '%% DEFAULT PARAMETER FILE\n\n');

% *** Default color scheme ***
% Note: I tried to make it somewhat colorbind friendly, but these might need to be
% tweaked for different users
% Colorscheme selected from https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=8
colorMat = [ 153, 153, 153;     % acoustic data (grey)
    0, 0, 0;            % Unlabeled detections (black)
    251, 154, 153;      % Whale 1 (pink)
    227, 26, 28;        % Whale 2 (red)
    193, 223, 138;      % Whale 3 (light green)
    56, 185, 7;         % Whale 4 (dark green)
    253, 191, 111;      % Whale 5 (orange cream)
    31, 120, 180;       % Whale 6 (blue)
    166, 206, 227;      % Whale 7 (light blue)
    106, 61, 154;       % Whale 8 (purple)
    202, 178, 214];     % Buzzes (lavender)
colorMat = colorMat./255; % convert from RGB to Matlab color scale
brushing.params.colorMat = colorMat;

fprintf(fid, '%% Color scheme:\n');
fprintf(fid, 'brushing.params.colorMat(1, :) = [%f, %f, %f]; %% acoustic data \n', colorMat(1, :));
fprintf(fid, 'brushing.params.colorMat(2, :) = [%f, %f, %f]; %% Unlabeled detections \n', colorMat(2, :));
for icol = 3:length(colorMat)-1
    fprintf(fid, 'brushing.params.colorMat(%d, :) = [%f, %f, %f]; %% Whale %d \n', icol, colorMat(icol, :), icol-2);
end

fprintf(fid, 'brushing.params.colorMat(%d, :) = [%f, %f, %f]; %% Buzzes \n', length(colorMat), colorMat(end, :));

% *** Figure Positions *** 
% Legends:
brushing.params.colorLegendPos = [1, 30, 175, 500]; % Label Color Legend figure position
fprintf(fid, '\nbrushing.params.colorLegendPos = [1, 30, 200, 500]; %% Label Color Legend Position \n');
brushing.params.commandLegendPos = [1, 530, 175, 300]; % keyboard command Legend figure position
fprintf(fid, 'brushing.params.commandLegendPos = [1, 560, 175, 300]; %% Keyboard command Legend Position \n');

% Main figure:
% calculate default figure size
fig = uifigure('name', 'Generating test figure'); % generate random figure
fig.WindowState = 'maximized'; % maximize the figure
drawnow
pause(.5) % MATLAB likes to execute the next line before the figure has finished maximizing the figure, so you have to add a pause so it has time to finish
maxFigPosition = fig.Position; % get position of maximized figure
close(fig)

xPosIn = brushing.params.colorLegendPos(3)+brushing.params.colorLegendPos(1); % find position of right side of legend

brushing.params.figPosition = maxFigPosition;
brushing.params.figPosition(1) = xPosIn; % begin brushDet figure to the right of legend
brushing.params.figPosition(3) = maxFigPosition(3)-xPosIn; % Correct figure width so right side of brushDet ends at edge of screen

fprintf(fid, '\nbrushing.params.figPosition = [%d, %d, %d, %d]; %% Figure position\n', brushing.params.figPosition);

fclose(fid);

end
