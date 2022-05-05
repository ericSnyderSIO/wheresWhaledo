function makeParamFile
% make a parameter file for function settings

global brushing

fid = fopen('brushing.params', 'w');
fprintf(fid, '%% DEFAULT PARAMETER FILE\n\n');
fprintf(fid, 'global brushing %% initialize global variable\n\n'); 

% *** Default color scheme ***
% Note: I tried to make it somewhat colorbind friendly, but these might need to be
% tweaked for different users
% Colorscheme selected from https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=8
colorMat = [ 153, 153, 153;     % acoustic data (grey)
    0, 0, 0;            % Unlabeled detections (black)
    251, 154, 153;      % Whale 1 (pink)
    193, 223, 138;      % Whale 2 (light green)
    166, 206, 227;      % Whale 3 (light blue)
    253, 191, 111;      % Whale 4 (orange cream)
    31, 120, 180;       % Whale 5 (blue)
    106, 61, 154;       % Whale 6 (purple)
    56, 185, 7;         % Whale 7 (dark green)
    227, 26, 28;        % Whale 8 (red)
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
brushing.params.figPosition(4) = brushing.params.figPosition(4) - 57;

fprintf(fid, '\nbrushing.params.figPosition = [%d, %d, %d, %d]; %% Figure position\n', brushing.params.figPosition);

fclose(fid);

end