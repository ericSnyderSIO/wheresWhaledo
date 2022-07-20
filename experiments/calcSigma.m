% Hydrophone location data (Lat, Lon, Depth, RMS)

% EE:
EEdep = [32.65871	-119.47711	-1325.5285	4.1544];    % deployment location
EErec = [32.65879	-119.47705	-1319.6305	3.1279];    % recovery location
[EEx, EEy] = latlon2xy_wgs84(EEdep(1), EEdep(2), EErec(1), EErec(2)); % difference between deployment and recover
EEdiff = sqrt(EEx^2 + EEy^2 + (EEdep(3)-EErec(3)).^2); % distance between dep and rec
EEloc = mean([EEdep(1:3); EErec(1:3)]); % Location used in model
sigEE_h = sqrt(EEdep(4)^2 + EErec(4)^2 + (EEdiff/2)^2); % std of hydrophone location

% EW:
EWdep = [32.65646	-119.48815	-1330.1631	14.4774];
EWrec = [32.65632	-119.48814	-1323.6842	3.8791];
[EWx, EWy] = latlon2xy_wgs84(EWdep(1), EWdep(2), EWrec(1), EWrec(2));
EWdiff = sqrt(EWx^2 + EWy^2 + (EWdep(3)-EWrec(3)).^2);
EWloc = mean([EWdep(1:3); EWrec(1:3)]);
sigEW_h = sqrt(EWdep(4)^2 + EWrec(4)^2 + (EWdiff/2)^2); % std of hydrophone location

% EN:
ENdep = [32.66219	-119.48425	-1333.0526	4.4839];
ENrec = [32.66221	-119.48424	-1321.2775	4.0117];
[ENx, ENy] = latlon2xy_wgs84(ENdep(1), ENdep(2), ENrec(1), ENrec(2));
ENdiff = sqrt(ENx^2 + ENy^2 + (ENdep(3)-ENrec(3)).^2);
ENloc = mean([ENdep(1:3); ENrec(1:3)]);
sigEN_h = sqrt(ENdep(4)^2 + ENrec(4)^2 + (ENdiff/2)^2); % std of hydrophone location

% ES:
ESdep = [32.65345	-119.48455	-1328.9836	13.5117];
ESrec = [32.65352	-119.48446	-1331.3959	6.53];
[ESx, ESy] = latlon2xy_wgs84(ESdep(1), ESdep(2), ESrec(1), ESrec(2));
ESdiff = sqrt(ESx^2 + ESy^2 + (ESdep(3)-ESrec(3)).^2);
ESloc = mean([ESdep(1:3); ESrec(1:3)]);
sigES_h = sqrt(ESdep(4)^2 + ESrec(4)^2 + (ESdiff/2)^2); % std of hydrophone location

[h(1, 1), h(1, 2)] = latlon2xy_wgs84(EWloc(1), EWloc(2), EEloc(1), EEloc(2)); % EE location in meters
h(1, 3) = EWloc(3)-EEloc(3);

[h(2, 1), h(2, 2)] = latlon2xy_wgs84(ENloc(1), ENloc(2), EEloc(1), EEloc(2)); % EE location in meters
h(2, 3) = ENloc(3)-EEloc(3);

[h(3, 1), h(3, 2)] = latlon2xy_wgs84(ESloc(1), ESloc(2), EEloc(1), EEloc(2)); % EE location in meters
h(3, 3) = EWloc(3)-ESloc(3);

note = 'EE is (0,0,0). Each column is [x,y,z] of a HARP. Rows are 1. EW; 2. EN; 3. ES.';
save('D:\SOCAL_E_63\xwavTables\instrumentLocs.mat', 'h', 'note') 



%% Sound Speed at depth
% calculated in D:\SOCAL_E_63\tracking\experiments\soundSpeedADCP.m

c = 1488.4; % sound speed 
sig_csml = 0.1345; % std of sound speed

%% Large Ap travel time due to uncertainty in sound speed
% calculated in
% D:\MATLAB_addons\gitHub\wheresWhaledo\experiments\rayBendingError

sig_clrg = .0066;

%% Small app TDOA uncertainty due to ray bending
% calculated in
% D:\MATLAB_addons\gitHub\wheresWhaledo\experiments\rayBendingError
sigEE_rayBend = 1.6858e-05;
sigEW_rayBend = 1.6111e-05;


%% Small app hydrophone position uncertainty
% calculated in:
% D:\SOCAL_E_63\tracking\experiments\inverseShipLocalization\arrayRPYfromGPS_usingTDOAs_experimentingWithInversion_EE_210518
% D:\SOCAL_E_63\tracking\experiments\inverseShipLocalization\arrayRPYfromGPS_usingTDOAs_experimentingWithInversion_EW_210518

sigEE_Hx = 0.1753;
sigEE_Hy = 0.1696;
sigEE_Hz = 0.0697;
sigEE_H = sqrt(sigEE_Hx^2 + sigEE_Hy^2 + sigEE_Hz^2);

sigEW_Hx = 0.1838;
sigEW_Hy = 0.1855;
sigEW_Hz = 0.1236;
sigEW_H = sqrt(sigEW_Hx^2 + sigEW_Hy^2 + sigEW_Hz^2);

%% Uncertainty due to drift
% Calculated in
% D:\SOCAL_E_63\tracking\experiments\clockSync\comeUpWithDriftForFullDeployment

load('D:\SOCAL_E_63\xwavTables\drift')

sigEW_drift = drift_std(1);
sigEN_drift = drift_std(2);
sigES_drift = drift_std(3);
%% 
% sigEE_equation = 'sigEE = sqrt(sigEE_H^2 + (TDOA*sig_csml)^2 + c^2*(sigEE_xcorr + sigEE_rayBend)^2)';
% sigEW_equation = 'sigEW = sqrt(sigEW_H^2 + (TDOA*sig_csml)^2 + c^2*(sigEW_xcorr + sigEW_rayBend)^2)';
% sigLargeAp_equation = 'sig = sqrt(sigEX_h^2 + sigEY_h^2 + (TDOA+drift)^2sig_clrg + c^2*(sigEX_EY_xcorr^2 + sigEX_drift^2+sigEY_drift^2))';

sigEE_equation = 'sigEE = sqrt(sig2EE*[ones(1,6); TDOA.^2; ones(1,6); sig_xcov.^2])';
sig2EE(1) = sigEE_H^2;      
sig2EE(2) = sig_csml^2;
sig2EE(3) = c^2*sigEE_rayBend^2;
sig2EE(4) = c^2; 

sigEW_equation = 'sigEW = sqrt(sig2EW*[ones(1,6); TDOA.^2; ones(1,6); sig_xcov.^2])';
sig2EW(1) = sigEW_H^2;      
sig2EW(2) = sig_csml^2;
sig2EW(3) = c^2*sigEE_rayBend^2;
sig2EW(4) = c^2; 

% large aperture variance:
sigLargeAp_equation = 'sig_lrg = sqrt(sig2lrg*[ones(1,6), (TDOA+drift).^2, ones(1,6), sig_xcov^2])';
sigLarge_notes = 'TDOA must be 1'
sig2lrg(1, 1) = sigEE_h^2 + sigEW_h^2;
sig2lrg(1, 2) = sigEE_h^2 + sigEN_h^2;
sig2lrg(1, 3) = sigEE_h^2 + sigES_h^2;
sig2lrg(1, 4) = sigEW_h^2 + sigEN_h^2;
sig2lrg(1, 5) = sigEW_h^2 + sigES_h^2;
sig2lrg(1, 6) = sigEN_h^2 + sigES_h^2;

sig2lrg(2, 1:6) = sig_clrg^2;

sig2lrg(3, 1) = c^2*sigEW_drift^2;
sig2lrg(3, 2) = c^2*sigEN_drift^2;
sig2lrg(3, 3) = c^2*sigES_drift^2;
sig2lrg(3, 4) = c^2*(sigEW_drift^2 + sigEN_drift^2);
sig2lrg(3, 5) = c^2*(sigEW_drift^2 + sigES_drift^2);
sig2lrg(3, 6) = c^2*(sigEN_drift^2 + sigES_drift^2);

sig2lrg(4, 1:6) = c^2;


save('sigmaValues', 'sigLargeAp_equation', 'sigEE_equation', 'sigEW_equation', 'sig2EE', 'sig2EW', 'sig2lrg')