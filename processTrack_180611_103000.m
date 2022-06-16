clear all
% Set up data and whatnot

encounterStart = datenum([18 06 11 10 30 00.1]);
encounterEnd = datenum([18 06 11 12 30 00]);

encName = '180611_1030';

% % load in DOA data:
% DOA{1} = load('D:\SOCAL_E_63\tracking\forRyan\SOCAL_E_63_EE\SOCAL_E_63_EE_C4_disk07_TDOA_DOA.mat');
% DOA{2} = load('D:\SOCAL_E_63\tracking\forRyan\SOCAL_E_63_EW\SOCAL_E_63_EW_C4_disk07_TDOA_DOA.mat');
% 
% for idoa = 1:numel(DOA)
%     Idoa{idoa} = find(DOA{idoa}.TDet>=encounterStart & DOA{idoa}.TDet<=encounterEnd);
% end

c = 1488.4;
fs1 = 200e3;
fs4 = 100e3;
spd = 60*60*24;

% filter parameters
fc = 20e3;
[b4, a4] = ellip(4,0.1,40,fc*2/fs4,'high'); % 4ch filter coeff's
[b1, a1] = ellip(4,0.1,40,fc*2/fs1,'high'); % 1ch filter coeff's

%% Detect clicks on all instruments

% load in xwav tables
XH{1} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EE_C4_xwavLookupTable.mat');
XH{2} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EW_C4_xwavLookupTable.mat');
XH{3} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EN_xwavLookupTable.mat');
XH{4} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_ES_xwavLookupTable.mat');

%% load in H matrices

hyd1 = load('D:\SOCAL_E_63\tracking\experiments\inverseProblem\matfiles\SOCAL_E_63_EE_Hmatrix_fromHydLocInversion_210702.mat');
hyd2 = load('D:\SOCAL_E_63\tracking\experiments\inverseProblem\matfiles\SOCAL_E_63_EW_Hmatrix_fromHydLocInversion_210702.mat');

% HEW = H;

% Reorder hydrophones to fit new TDOA order
H{1} = [hyd1.hydPos(2,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(3,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(3,:)-hyd1.hydPos(2,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(2,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(3,:)];
     
H{2} = [hyd2.hydPos(2,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(3,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(3,:)-hyd2.hydPos(2,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(2,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(3,:)];


%% Run Detector
% load('detections.mat')
% 
% SOCAL_E_63_EE
[DET{1}] = detectClicks_4ch_SOCAL_E_63(encounterStart, encounterEnd, XH{1}.xwavTable, H{1}, c, 'detClicks_4ch.params'); 

% SOCAL_E_63_EW
[DET{2}] = detectClicks_4ch_SOCAL_E_63(encounterStart, encounterEnd, XH{2}.xwavTable, H{2}, c, 'detClicks_4ch.params');
 
% SOCAL_E_63_EN
[DET{3}] = detectClicks_1ch(encounterStart, encounterEnd, XH{3}.xwavTable, 'detClicks_1ch.params');

% SOCAL_E_63_ES
[DET{4}] = detectClicks_1ch(encounterStart, encounterEnd, XH{4}.xwavTable, 'detClicks_1ch.params');
 
save(['detections', encName], 'DET')


%% run brushDOA
load('detections.mat')

[det1, det2] = brushDOA(DET{1}, DET{2});
% 
DET{1} = det1;
DET{2} = det2;
 
save(['detections_brushDOA', encName, '.mat'], 'DET')

%% Click train correlation for whole encounter, using whale labels

