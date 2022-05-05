function [T, TDOA] = tdoaFromShipSound(tstart, tend, xwavTable, varargin)
% [T, TDOA] = tdoaFromShipSound(tstart, tend, xwavTable)
% [T, TDOA] = tdoaFromShipSound(tstart, tend, xwavTable, paramsFile)
% calculates the Time Difference of Arrival by cross-covariating the ship
% sound.
% VARIABLES:
% - tstart: beginning of ship passage (in days since 01-Jan-2000 00:00:00)
% - tend: end of ship passage (in days since 01-Jan-2000 00:00:00)
% - xwavTable: the table of xwav start times and paths to files
%       can be generated using the readxwavSegment remora for Triton
% - OPTIONAL: paramsFile is the path and name of the .params file containing
%       all necessary parameters. Without this, user will be prompted to
%       select the .params file


if nargin==3 % No parameter file included: prompt user to select one
    [file, path] = uigetfile('*.params', 'Select a .params file for tdoaFromShipSound');
    
    loadParams(fullfile(path, file))


elseif nargin==4
    loadParams(varargin{1})
end

global shipTDOA % global variable that pulls in params 

% determine filter coefficients
[b, a] = ellip(4,0.1,40,shipTDOA.fc.*2/shipTDOA.fs);

Nxc = shipTDOA.xcovSegmentSize * shipTDOA.fs; % number of samples to use in xcov

spd = 60*60*24; % number of seconds per day
Ntdoa = round((tend-tstart)*spd/shipTDOA.xcovSegmentSize); % Estimate number of TDOAs that will be computed

% initialize variables:
T = zeros(Ntdoa, 1);
TDOA = zeros(Ntdoa, length(shipTDOA.ixcov));

% load time period:
tLoad(1) = tstart; 
tLoad(2) = tstart + shipTDOA.loadSegmentSize/spd;

i = 0; % counter

while tLoad(2) <= tend
    fprintf('\nCalculating TDOA for segment %d of %d\n', i+1, Ntdoa)
    [x, t] = quickxwavRead(tLoad(1), tLoad(2), shipTDOA.fs, xwavTable); % read in xwav data
    xf = filtfilt(b, a, x); % filter data
    
    % samples to use in xcov:
    nxc(1) = 1; 
    nxc(2) = nxc(1) + Nxc;

    while nxc(2) <= length(xf)
        i = i+1;
        T(i) = (t(nxc(2)) + t(nxc(1)))/2; % time stamp is center of xcov segment
        
        xseg = xf(nxc(1):nxc(2), :);

        [xc, lags] = xcov(xseg);

        for npair = 1:length(shipTDOA.ixcov)
            [~, m] = max(xc(:, shipTDOA.ixcov(npair)));
            TDOA(i, npair) = lags(m)/shipTDOA.fs;
        end

        nxc(1) = nxc(2) + 1;
        nxc(2) = nxc(1) + Nxc;
    end
    
    tLoad(1) = tLoad(2);       
    tLoad(2) = tLoad(1) + shipTDOA.loadSegmentSize/spd;
end


% remove excess zeros from initialization

Irem = find(T==0);

T(Irem) = [];
TDOA(Irem, :) = [];