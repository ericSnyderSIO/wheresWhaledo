clear all
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

h0 = mean([EEloc; EWloc]);

[h(1, 1), h(1, 2)] = latlon2xy_wgs84(EEloc(1), EEloc(2), h0(1), h0(2)); % EE location in meters
h(1, 3) = EEloc(3)-h0(3);

[h(2, 1), h(2, 2)] = latlon2xy_wgs84(EWloc(1), EWloc(2), h0(1), h0(2)); % EE location in meters
h(2, 3) = EWloc(3)-h0(3);

[h(3, 1), h(3, 2)] = latlon2xy_wgs84(ENloc(1), ENloc(2), h0(1), h0(2)); % EE location in meters
h(3, 3) = ENloc(3)-h0(3);

[h(4, 1), h(4, 2)] = latlon2xy_wgs84(ESloc(1), ESloc(2), h0(1), h0(2)); % EE location in meters
h(4, 3) = EWloc(3)-h0(3);

hLatLonZ = [EEloc(1:3); EWloc(1:3); ENloc(1:3); ESloc(1:3)]

hyd1a = load('D:\SOCAL_E_63\tracking\experiments\inverseProblem\matfiles\SOCAL_E_63_EE_Hmatrix_fromHydLocInversion_210702.mat');
hyd1b = load('D:\MATLAB_addons\gitHub\wheresWhaledo\receiverPositionInversion\SOCAL_E_63_EE_Hmatrix_new.mat');
% hyd1.hydPos = hyd1.recPos;
hyd2a = load('D:\SOCAL_E_63\tracking\experiments\inverseProblem\matfiles\SOCAL_E_63_EW_Hmatrix_fromHydLocInversion_210702.mat');
hyd2b = load('D:\MATLAB_addons\gitHub\wheresWhaledo\receiverPositionInversion\SOCAL_E_63_EW_Hmatrix_new.mat');

% Reorder hydrophones to fit new TDOA order
H1 = [hyd1a.hydPos(2,:)-hyd1a.hydPos(1,:);
    hyd1a.hydPos(3,:)-hyd1a.hydPos(1,:);
    hyd1a.hydPos(4,:)-hyd1a.hydPos(1,:);
    hyd1a.hydPos(3,:)-hyd1a.hydPos(2,:);
    hyd1a.hydPos(4,:)-hyd1a.hydPos(2,:);
    hyd1a.hydPos(4,:)-hyd1a.hydPos(3,:)];

H2 = [hyd2a.hydPos(2,:)-hyd2a.hydPos(1,:);
    hyd2a.hydPos(3,:)-hyd2a.hydPos(1,:);
    hyd2a.hydPos(4,:)-hyd2a.hydPos(1,:);
    hyd2a.hydPos(3,:)-hyd2a.hydPos(2,:);
    hyd2a.hydPos(4,:)-hyd2a.hydPos(2,:);
    hyd2a.hydPos(4,:)-hyd2a.hydPos(3,:)];

hloc = h;

% load drift:
load('D:\SOCAL_E_63\tracking\experiments\clockSync\drift.mat');
dp{1} = coeffvalues(Dpoly{1}); % drift coefficients between inst 1 and 2
dp{2} = coeffvalues(Dpoly{2}); % drift coefficients between inst 1 and 3
dp{3} = coeffvalues(Dpoly{3}); % drift coefficients between inst 1 and 4
dp{4} = coeffvalues(Dpoly{4}); % drift coefficients between inst 2 and 3
dp{5} = coeffvalues(Dpoly{5}); % drift coefficients between inst 2 and 4
dp{6} = coeffvalues(Dpoly{6}); % drift coefficients between inst 3 and 4

load('D:\SOCAL_E_63\tracking\interns2022\ericEdits_allTracks\track600_180611_110414\SOCAL_E_63_track600_180611_110414_ericMod_localized_cleaned.mat')
spd = 60*60*24;
for nd = 1:6
    drft(nd) = polyval(dp{nd}, mean(whale{1}.TDet));
end
M{1}.c = 1488.4;

%%
global LOC
loadParams('localize.params')

sig2sml = LOC.sig_sml^2; % variance of small ap
sig2lrg = LOC.sig_lrg^2; % variance of large ap
Asml = (2*pi*sig2sml)^(-12/2); % coefficient of small ap
Alrg = (2*pi*sig2lrg)^(-6/2); % coefficient of small ap

Iuse = find(sum(~isnan(whale{1}.TDOA), 2)>=7);

tstart = min(whale{1}.TDet(Iuse));
tend = max(whale{1}.TDet(Iuse));

ti = tstart:60/spd:tend;
iwhale{1}.TDet = ti;

for itdoa = 1:18
    Iuse = find(~isnan(whale{1}.TDOA(:, itdoa)));
    tdoa = smoothdata(whale{1}.TDOA(Iuse, itdoa));
    iwhale{1}.TDOA(:, itdoa) = interp1(whale{1}.TDet(Iuse), tdoa, ti, 'spline');
end

% TDOAi(:, 13:end) = TDOAi(:, 13:end) - drft;

Asml = (2*pi*sig2sml)^(-12/2); % coefficient of small ap
Alrg = (2*pi*sig2lrg)^(-6/2); % coefficient of large ap

for i = 1:length(ti)

    whaleOut = localize(iwhale, hloc, H1, H2, dp);

end

%%

% Irem = find(wi(:, 2)<-2000);
% wi(Irem, :) = [];
% ti(Irem) = [];
% TDOAi(Irem, :) = [];

figure(12)
plot3(whaleOut{1}.wloc(:, 1), whaleOut{1}.wloc(:, 2), whaleOut{1}.wloc(:, 3), '-')
hold on
scatter3(whaleOut{1}.wloc(:, 1), whaleOut{1}.wloc(:, 2), whaleOut{1}.wloc(:, 3), [], (ti-min(ti))./max(ti))
hold off


figure(13)
for sp = 1:3
    subplot(3,1,sp)
    plot(ti, whaleOut{1}.wloc(:, sp))
    hold on
    scatter(ti, whaleOut{1}.wloc(:, sp), [], (ti-min(ti))./max(ti))
    hold off
    datetick
end

figure(14)
for sp = 1:18
    subplot(3,6,sp)
    plot(iwhale{1}.TDOA(:,sp), '.')
end