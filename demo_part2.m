% close all force
% clear all
clear inDet DATA


%% ******************** DEMO, PART 2: Real data and subsequent functions ************************

tCurrent = datenum([2018, 7, 1, 19, 14, 22]);

% xwav tables:
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EE_C4_xwavLookupTable');
XH{1} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EW_C4_xwavLookupTable');
XH{2} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EN_xwavLookupTable');
XH{3} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_ES_xwavLookupTable');
XH{4} = xwavTable;

% 4 channel parameters
H{1} = load('D:\SOCAL_E_63\tracking\experiments\largeApertureTDOA\wheresWhaledo2.0\receiverPositionInversion\SOCAL_E_63_EE_Hmatrix.mat');

% *************CHANGE TO EW WHEN YOU HAVE IT!!!!!!!!!! *************
H{2} = load('D:\SOCAL_E_63\tracking\experiments\largeApertureTDOA\wheresWhaledo2.0\receiverPositionInversion\SOCAL_E_63_EE_Hmatrix.mat');


% various parameters:
spd = 60*60*24; % seconds per day, for converting days from datenum to seconds
y2k = datenum([2000, 0, 0, 0, 0, 0]); % for changing 2018 to 0018 and vice versa
c = 1488.4;

fs4 = 100e3; % 4ch sampling rate
fs1 = 200e3; % single channel sampling rate

% drift correction
load('D:\SOCAL_E_63\xwavTables\drift') % drift (sec) and tdrift (datenum) for EW, EN, ES relative to EE
tdrifti = tdrift(1):60/spd:tdrift(end);
for nd = 1:3
    drifti(:,nd) = interp1(tdrift, drift(nd, :).', tdrifti);
end

% filter parameters
fc = 20e3;
[b4, a4] = ellip(4,0.1,40,fc*2/fs4,'high'); % 4ch filter coeff's
[b1, a1] = ellip(4,0.1,40,fc*2/fs1,'high'); % 1ch filter coeff's

%% step 1) load in 30s of data from all four HARPs, highpass filter

% will eventually need this in a while loop with a few keyboard command
% options (maybe dialogue box?)

% start at beginning of track and work towards the earlier detections
clear global
tstart = tCurrent(1)-y2k-45/spd;
tend = tstart + 30/spd;
TEND = tend + 55/spd;

% while tend<TEND
clear X inDet
[xee, tee] = quickxwavRead(tstart, tend, fs4, XH{1});
xee = filtfilt(b4, a4, xee);
xee = xee.';

[xew, tew] = quickxwavRead(tstart, tend, fs4, XH{2});
id = find(abs(tdrifti-mean(tew))==min(abs(tdrifti-mean(tew))));
meanDrift(1) = drifti(id, 1);
tew = tew - drifti(id, 1)/spd; % correct for drift
xew = filtfilt(b4, a4, xew);
xew = xew.';

[xen, ten] = quickxwavRead(tstart, tend, fs1, XH{3});
id = find(abs(tdrifti-mean(ten))==min(abs(tdrifti-mean(ten))));
meanDrift(2) = drifti(id, 2);
ten = ten - drifti(id, 2)/spd; % correct for drift
xen = filtfilt(b1, a1, xen);
xen = xen.';

[xes, tes] = quickxwavRead(tstart, tend, fs1, XH{4});
id = find(abs(tdrifti-mean(tes))==min(abs(tdrifti-mean(tes))));
meanDrift(3) = drifti(id, 3);
tes = tes - drifti(id, 3)/spd; % correct for drift
xes = filtfilt(b1, a1, xes);
xes = xes.';

% put data into format for brushDet
X{1} = [(tee-tstart).*spd; xee(1,:)];
X{2} = [(tew-tstart).*spd; xew(1,:)];
X{3} = [(ten-tstart).*spd; xen];
X{4} = [(tes-tstart).*spd; xes];

% put data into format for brushDOA
X4ch{1} = [(tee-tstart).*spd; xee];
X4ch{2} = [(tew-tstart).*spd; xew];

rowxcov = [2, 3, 4, 7, 8, 12]; % which rows of the xcov do we need for TDOA

% detect clicks

thresholds = [25,25, 35, 35];
for irec = 1:numel(X)
    [~, idet{irec}] = findpeaks(X{irec}(2,:), 'MinPeakDistance', 800, 'MinPeakHeight', thresholds(irec));

    inDet{irec}(1,:) = X{irec}(1, idet{irec});
    inDet{irec}(2,:) = X{irec}(2, idet{irec});
end

[outDet, labels] = brushDet_old(X, inDet)

%% DOA clean up code
% for SOCAL_E_63_EE
% temp = inputdlg('Which array had the best detections?');
% primaryArray = str2double(temp{1});

Hmat{1} = H{1}.H;
Hmat{2} = H{2}.H;

desLabel = '0'; % label to examine
Nwin = 200; % window over which to xcov

[TDet, DOA, TDOA, newlabels] = brushDOA(X4ch, outDet, labels, Hmat, desLabel, Nwin, rowxcov, fs4, c);



% Click train correlation


% tstart = tend;
% tend = tstart + 30/spd;
% end