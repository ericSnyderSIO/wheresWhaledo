% An example of how to invert for hydrophone positions in a small-aperture 
% array from the TDOA of known source locations (ship engine sounds)
close all
%%
deploymentName = 'SOCAL_E_63_EE';

XH = load(['D:\SOCAL_E_63\xwavTables\', deploymentName, '_C4_xwavLookupTable']);

ymd = [18 03 16]; % year, month, day of ship localization period (should be in metadata)

filename = [deploymentName, '_ttgps.txt'];
[Lat, Lon, T] = readttgps(filename, ymd);

save([deploymentName, '_shipLatLon'], 'Lat', 'Lon', 'T')

%% Calculate TDOA
% load([deploymentName, '_shipLatLon.mat'])
tstart = T(1);
tend = T(end);

% [T, TDOA] = tdoaFromShipSound(tstart, tend, XH.xwavTable, 'tdoaFromShipSound.params');
% [T, TDOA] = tdoaFromShipSound(tstart, tend, XH.xwavTable);


% save([deploymentName, '_shipTDOA.mat'], 'T', 'TDOA')

%% Brush TDOA and invert for hydrophone positions
TDOAstruct = load([deploymentName, '_shipTDOA.mat']);
shipStruct = load([deploymentName, '_shipLatLon.mat']);

% close all
% H = brushDet_hydPosInv('shipTDOA.mat')
[H, hydPos, CI] = brushTDOA_hydPosInv(TDOAstruct, shipStruct, 'brushTDOA_hydPosInv_SOCAL_E_63_EE.txt');




