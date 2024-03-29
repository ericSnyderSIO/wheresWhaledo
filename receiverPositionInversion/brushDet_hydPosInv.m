function [H, recPos, CI] = brushDet_hydPosInv(varargin)
% Runs GUI to brush TDOA detections. Calculates the hydrophone positions
% and H matrix.
% H = brushDet_hydPosInv prompts user to select files for TDOA, ship
% positions, and parameter file.
% 
% brushDet_hydPosInv('TDOA.mat', 'shipLoc.mat') loads data from TDOA.mat
% and shipLoc.mat files and prompts user to select parameter file.
%
% brushDet_hydPosInv(TDOAstruct, shipLocStruct) where TDOA and shipLoc are 
% MATLAB structs assumes TDOAstruct.TDOA is NdetxNrec matrix of TDOAs, and
% TDOAstruct.T is the associated time vector. 
% shipLocStruct.Lat, shipLocStruct.Lon, and shipLocStruct.T must each be 
% Nx1 vectors of the Latitude, Longitude, and Time of the ship locations.
%
% brushDet_hydPosInv(..., paramFile) loads in the parameters in the file 
% specified by paramFile.
%
% OUTPUTS: H is the H matrix of vectors drawn between each receiver pair
% recPos are the receiver positions, where receiver 1 is (0,0,0)
% CI are the 95% confidence intervals obtained by inverting for recPos.

%% Load in data
switch nargin
    case 0
        % select TDOA file
        [file, path] = uigetfile('.mat', 'Select TDOA file', 'Select TDOA file');
        TDOAstruct = load(fullfile(path, file));

        TDOA = TDOAstruct.TDOA; % TDOA of acoustic data
        tAc = TDOAstruct.T; % time of acoustic data

        % select ship position file
        [file, path] = uigetfile('.mat', 'Select ship positions file', 'Select ship positions file');

        shipLocStruct =  load(fullfile(path, file));
        Lat = shipLocStruct.Lat;
        Lon = shipLocStruct.Lon;
        tShip = shipLocStruct.t; % time of ship data

        % select parameter file
        [file, path] = uigetfile('.mat', 'Select ship positions file', 'Select ship positions file');
        loadParams(fullfile(path, file))

    case 1 % Time and TDOA included; no ship positions or parameter file

        if isa(varargin{1}, 'char') % user provided path to TDOA file
            TDOAstruct = load(varargin{1});

            TDOA = TDOAstruct.TDOA; % TDOA of acoustic data
            tAc = TDOAstruct.T; % time of acoustic data

        elseif isa(varargin{1}, 'struct')
            TDOA = varargin{1}.TDOA; % TDOA of acoustic data
            tAc = varargin{1}.T; % time of acoustic data

        else
            errordlg('Invalid input type for TDOA', 'Error')

        end

        % select ship position file
        [file, path] = uigetfile('.mat', 'Select ship positions file', 'Select ship positions file');

        shipLocStruct =  load(fullfile(path, file));
        Lat = shipLocStruct.Lat;
        Lon = shipLocStruct.Lon;
        tShip = shipLocStruct.t; % time of ship data

        % select parameter file
        [file, path] = uigetfile('.mat', 'Select ship positions file', 'Select ship positions file');
        loadParams(fullfile(path, file))

    case 2 % Time and TDOA included; ship positions included; no parameter file

        if isa(varargin{1}, 'char') % user provided path to TDOA file
            TDOAstruct = load(varargin{1});

            TDOA = TDOAstruct.TDOA; % TDOA of acoustic data
            tAc = TDOAstruct.T; % time of acoustic data

        elseif isa(varargin{1}, 'struct') % user provided TDOA struct
            TDOA = varargin{1}.TDOA; % TDOA of acoustic data
            tAc = varargin{1}.T; % time of acoustic data

        else
            errordlg('Invalid input type for TDOA', 'Error')

        end

        if isa(varargin{2}, 'char') % user provided path to ship lat/lon/t file
            shipLocStruct = load(varargin{2});

            Lat = shipLocStruct.Lat;
            Lon = shipLocStruct.Lon;
            tShip = shipLocStruct.T; % time of ship data

        elseif isa(varargin{2}, 'struct') % user provided  to ship lat/lon/t struct
            Lat = varargin{2}.Lat;
            Lon = varargin{2}.Lon;
            tShip = varargin{2}.T; % time of ship data

        else
            errordlg('Invalid input type for ship lat/lon/z', 'Error')

        end


        % select parameter file
        [file, path] = uigetfile('.params', 'Select parameter file');
        loadParams(fullfile(path, file))

    case 3 % all variables included
        if isa(varargin{1}, 'char') % user provided path to TDOA file
            TDOAstruct = load(varargin{1});

            TDOA = TDOAstruct.TDOA; % TDOA of acoustic data
            tAc = TDOAstruct.T; % time of acoustic data

        elseif isa(varargin{1}, 'struct') % user provided TDOA struct
            TDOA = varargin{1}.TDOA; % TDOA of acoustic data
            tAc = varargin{1}.T; % time of acoustic data

        else
            errordlg('Invalid input type for TDOA', 'Error')

        end

        if isa(varargin{2}, 'char') % user provided path to ship lat/lon/t file
            shipLocStruct = load(varargin{1});

            Lat = shipLocStruct.Lat;
            Lon = shipLocStruct.Lon;
            tShip = shipLocStruct.T; % time of ship data

        elseif isa(varargin{2}, 'struct') % user provided  to ship lat/lon/t struct
            Lat = varargin{2}.Lat;
            Lon = varargin{2}.Lon;
            tShip = varargin{2}.T; % time of ship data

        else
            errordlg('Invalid input type for ship lat/lon/z', 'Error')

        end
       
        loadParams(varargin{3})

end

%% Set up necessary variables
Nrec = size(TDOA, 2); % number of receivers
[x, y] = latlon2xy(Lat, Lon, )
% time frame for analysis
t1 = max([tShip(1), tAc(1)]); 
t2 = min([tShip(end), tAc(end)]); 

spd = 60*60*24; % seconds per day
ti = t1:1/spd:t2; % set time vector for interpolation

% Solving 

%% Run GUI
str = 'c';

while ~strcmp(str, 'q')
    
end

[file,path] = uiputfile('Select location to save hydrophone parameters')

save(fullfile(path, file), 'H', 'rLatLonZ')