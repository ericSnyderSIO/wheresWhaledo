function [TDet, DOA, TDOA, labels] = brushDOA(X, det, labels, H, desLabel, Nwin, rowxcov, fs, c, varargin)
% runs GUI to calculate direction of Arrival in 4ch data
% X{narr}(1,:) is time vector for array 'narr'
% X{narr}(nch+1,:) is acoustic data on channel 'nch' of array 'narr'
% det{narr}(1, :) are the detection times on array 'narr'
% det{narr}(2, :) are the detection amplitudes on arra 'narr'
% labels{narr} are the labels of each detection in det
% H{narr} is the H matrix of array 'narr'
% desLabel is the desired label string (i.e. '0' for unlabeled detections)
% Nwin is the window size around each detection used in the
% cross-covariance
% rowxcov are the rows of xcov to use to in the TDOA calculations
% fs is the sampling rate
% c is the speed of sound

global brushing

if nargin < 9
    errordlg('Not enough input arguments', 'brushDet Error')
elseif nargin == 9 % no parameter file specified
    % check for existing default param file
    if isfile('brushing.params')
        loadParams('brushing.params')
    else
        makeParamFile
    end
elseif nargin == 10 % parameter file specified
    loadParams(varargin{1})
else
    errordlg('incorrect number of input arguments', 'brushDet Error')
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

[TDet, DOA, TDOA, labels] = run_brushDOA(X, det, labels, H, desLabel, Nwin, rowxcov, fs, c);

qselect = input('\nEnter ''q'' to quit: ', 's');

end

function [TDet, DOA, TDOA, labels] = run_brushDOA(X, det, labels, H, desLabel, Nwin, rowxcov, fs, c)
global brushing

brushing.labels = labels;
brushing.Narr = numel(X); % number of arrays

% initialize figure
fig = findall(0, 'Type', 'figure', 'name', 'Brush DOA');
if isempty(fig)
    fig = figure('Name', 'Brush DOA');
end

fig.Position = brushing.params.figPosition;

% calculate DOA:
for iarr = 1:brushing.Narr
    brushing.Ilab{iarr} = find(labels{iarr}==desLabel); % Index of detections with desired label
    brushing.TDet{iarr} = det{iarr}(1, brushing.Ilab{iarr}); % Detection times

    [brushing.DOA{iarr}, brushing.TDOA{iarr}] = calcDOA(X{iarr}, det{iarr}(:, brushing. Ilab{iarr}), H{iarr}, Nwin, rowxcov, fs, c);
end

% plot DOA
brushing.Nsp = numel(X)*2;
isp = 0; % subplot number counter
for iarr = 1:numel(X)
    isp = isp+1;
    subplot(brushing.Nsp, 1, isp)
    scatter(brushing.TDet{iarr}, brushing.DOA{iarr}(1, :), 30, brushing.params.colorMat(brushing.colorIndex{iarr}, :), 'filled')
    title(['Array ', num2str(iarr), ' Azimuth'])
    ylabel('Angle (^\circ)')

    isp = isp+1;
    subplot(brushing.Nsp, 1, isp)
    scatter(brushing.TDet{iarr}, brushing.DOA{iarr}(2, :), 30, brushing.params.colorMat(brushing.colorIndex{iarr}, :), 'filled')
    title(['Array ', num2str(iarr), ' Elevation'])
    ylabel('Angle (^\circ)')
    xlabel('Time (s)')
end
brush on

hManager = uigetmodemanager(fig);
[hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2 (on 2014b or later)
set(fig, 'KeyPressFcn', @keyPressCallback);

DOA = brushing.DOA;
TDet = brushing.TDet;
labels = brushing.labels;
TDOA = brushing.TDOA;

end

%% Calculate DOA
function [DOA, TDOA] = calcDOA(X, det, H, Nwin, rowxcov, fs, c)
DOA = zeros(2, length(det));
TDOA = zeros(6, length(det));


for i = 1:length(det)
    tdet = det(1, i); % time of current detection
    [~, ind] = min(abs(X(1,:) - tdet)); % index where detection occurs in acoustic data

    i1 = max([1, ind-Nwin/2]); % beginning sample around detection
    i2 = min([ind+Nwin/2, length(X)]); % ending sample around detection

    x = X(2:5, i1:i2); % segment of acoustic data
    [xc, lags] = xcov(x.'); 

    tdoa = zeros([6, 1]); % initialize tdoa

    for itdoa = 1:length(rowxcov)
        [~, imax] = max(xc(:, rowxcov(itdoa)));
        tdoa(itdoa) = lags(imax)/fs;

    end

    TDOA(:, i) = tdoa.';
    s = -(tdoa.*c)\H;
    s = s./sqrt(sum(s.^2));

    DOA(1, i) = atan2d(s(2), s(1)); % azimuth
    DOA(2, i) = 180 - acosd(s(3)/sqrt(sum(s.^2))); % elevation
end

end

%% KeyPressCallback
function keyPressCallback(source, eventdata)
global brushing

key = eventdata.Key; % pressed key

% get highlighted points on each array
isp = 0; % subplot counter
for iarr = 1:brushing.Narr % iterate through arrays
    arrNum = brushing.Narr + 1 - iarr; % matlab stores subplots from bottom to top, so the array number will be flipped from iarr
    Ind{arrNum} = []; % Indices of selected points on array arrNo
    axNum{arrNum} = []; % axes numbers associated with array arrNo
    for idoa = 1:2 % iterate through subplots of each array
        isp = isp+1;
        selectedData = source.Children(isp).Children.BrushData;
        Itemp = find(selectedData~=0); % indices of selected data
        Ind{arrNum} = [Ind{arrNum}, Itemp];
        axNum{arrNum} = [axNum{arrNum}, isp];
    end

end

numkey = str2double(key); % convert keyboard input to a number (returns NaN if value is not a number)

if ~(isempty(numkey)||isnan(numkey))
    if numkey>8 || numkey<0
        errBox = msgbox('error: invalid whale number\nSelect a number 1 thorugh 8', 'Error');
    else % valid number selected
        for arrNum = 1:brushing.Narr
            brushing.labels{arrNum}(Ind{arrNum}) = key; % change labels
            brushing.colorIndex{arrNum}(Ind{arrNum}) = numkey+2;

            for iax = 1:length(axNum{arrNum})
                 
                set(source.Children(axNum{arrNum}(iax)).Children, 'cdata', brushing.params.colorMat(brushing.colorIndex{arrNum}(Ind{arrNum}), :) )
            end
        end
    end
else

end

end