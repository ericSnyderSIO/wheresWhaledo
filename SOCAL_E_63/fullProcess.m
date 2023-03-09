TSTART = datenum([18 06 11 10 20 00]);
TEND = datenum([18 06 11 12 20 00]);

c = 1488.4;

load('D:\SOCAL_E_63\xwavTables\instrumentLocs.mat')
hydLoc{1} = hLatLonZ(1,:);
hydLoc{2} = hLatLonZ(2,:);
hydLoc{3} = hLatLonZ(3,:);
hydLoc{4} = hLatLonZ(4,:);

h0 = mean([hydLoc{1}; hydLoc{2}]);

% convert hydrophone locations to meters:
[h1(1), h1(2)] = latlon2xy_wgs84(hydLoc{1}(1), hydLoc{1}(2), h0(1), h0(2));
h1(3) = abs(h0(3))-abs(hydLoc{1}(3));

[h2(1), h2(2)] = latlon2xy_wgs84(hydLoc{2}(1), hydLoc{2}(2), h0(1), h0(2));
h2(3) = abs(h0(3))-abs(hydLoc{2}(3));

[h3(1), h3(2)] = latlon2xy_wgs84(hydLoc{3}(1), hydLoc{3}(2), h0(1), h0(2));
h3(3) = abs(h0(3))-abs(hydLoc{3}(3));

[h4(1), h4(2)] = latlon2xy_wgs84(hydLoc{4}(1), hydLoc{4}(2), h0(1), h0(2));
h4(3) = abs(h0(3))-abs(hydLoc{4}(3));

hloc = [h1;h2;h3;h4];

% xwav tables:
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EE_C4_xwavLookupTable');
XH{1} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EW_C4_xwavLookupTable');
XH{2} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EN_xwavLookupTable');
XH{3} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_ES_xwavLookupTable');
XH{4} = xwavTable;

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

% load drift:
load('D:\SOCAL_E_63\tracking\experiments\clockSync\drift.mat');
dp{1} = coeffvalues(Dpoly{1}); % drift coefficients between inst 1 and 2
dp{2} = coeffvalues(Dpoly{2}); % drift coefficients between inst 1 and 3
dp{3} = coeffvalues(Dpoly{3}); % drift coefficients between inst 1 and 4
dp{4} = coeffvalues(Dpoly{4}); % drift coefficients between inst 2 and 3
dp{5} = coeffvalues(Dpoly{5}); % drift coefficients between inst 2 and 4
dp{6} = coeffvalues(Dpoly{6}); % drift coefficients between inst 3 and 4

fs(1:2) = 100e3;
fs(3:4) = 200e3;
%%
for n = 1:2
    DET{n} = detectClicks_4ch(TSTART, TEND, XH{n}, H{n}, c, 'D:\MATLAB_addons\gitHub\wheresWhaledo\detClicks_4ch.params')
end

for n = 3:4
    DET{n} = detectClicks_1ch(TSTART, TEND, XH{n}, 'D:\MATLAB_addons\gitHub\wheresWhaledo\detClicks_1ch.params')
end

%%

[DET{1}, DET{2}] = brushDOA(DET{1}, DET{2})

%%
colNum = length(unique(DET{1}.color))
for wn = 1:length(colNum)
    [CTC{wn}, DET]  = clickTrainCorr(DET, wn, [1,2], [1,2], [3,4]);
end
%% 

whale = calcTDOAfromCTC(CTC, XH)

%% Localize
whale = localize(whale, hloc, H{1}, H{2}, dp);
save('D:\SOCAL_E_63\tracking\interns2022\ericEdits_allTracks\track600_180611_110414\SOCAL_E_63_track600_180611_110414_ericMod_localized.mat', 'whale')
%% brush TDOA
load('D:\SOCAL_E_63\tracking\interns2022\ericEdits_allTracks\track600_180611_110414\SOCAL_E_63_track600_180611_110414_ericMod_localized.mat')
whaleOut = brushTDOA(whale, H);

%% smooth
clear sm
% Simplified:

LOCATIONS = whale{1}.wloc;
TIME = whale{1}.TDet;

Irem = find(isnan(LOCATIONS(:,1)));

LOCATIONS(Irem, :) = [];
TIME(Irem) = [];

Irem = find(TIME<6737.466 | LOCATIONS(:, 3)>1000 | LOCATIONS(:,2) <-2000 | LOCATIONS(:,1) <-700);

LOCATIONS(Irem, :) = [];
TIME(Irem) = [];

figure(98)
plot(TIME, LOCATIONS, '.')
datetick

sm(:,1) = smoothdata(LOCATIONS(:,1), 'movmean', 100);
sm(:,2) = smoothdata(LOCATIONS(:,2), 'movmean', 100);
sm(:,3) = smoothdata(LOCATIONS(:,3), 'movmean', 100);

figure(99)
plot(TIME, sm, '.'); datetick
hold on
plot(TIME, LOCATIONS, 'x')

LOCATIONS = sm;
save('oneWhaleExample', "TIME", 'LOCATIONS')
