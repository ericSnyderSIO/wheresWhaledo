c = 1488.4;
M = load('B:\TDOAmodel_50dx50dy15dz.mat');

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


fdir = dir('D:\SOCAL_E_63\tracking\interns2022\ericEdits_allTracks\track*');

ndet = 0;
ncluster = 0;
for ndf = 1:numel(fdir)
    fname = dir([fullfile(fdir(ndf).folder, fdir(ndf).name), '\*_localized.mat']);
    if isempty(fname)
        continue
    end
    load(fullfile(fname.folder, fname.name));

    for wn = 1:numel(whale)
        Iuse = find(sum(~isnan(whale{wn}.TDOA),2)==18);
        nstart = ndet+1;
        if ~isempty(Iuse)
            for i = 1:length(Iuse)
                % localize with DOAintersect
                tdoa1 = whale{wn}.TDOA(Iuse(i), 1:6);
                doa1 = (tdoa1.'\H{1})./c;
                doa1 = doa1./sqrt(sum(doa1.^2));

                tdoa2 = whale{wn}.TDOA(Iuse(i), 7:12);
                doa2 = (tdoa2.'\H{2})./c;
                doa2 = doa2./sqrt(sum(doa2.^2));

                D = [doa1; -doa2];
                R = D.'\(h2-h1).';

                w1 = R(1).*doa1 + h1;
                w2 = R(2).*doa2 + h2;

                wloc = mean([w1; w2]);

                if 1%abs(wloc(1))<2000 & abs(wloc(2))<3000 & abs(wloc(2))>300 & wloc(3)<800
                    ndet = ndet+1;
                    E.wloc_doa(ndet, :) = wloc;
                    E.wloc_mod(ndet, :) = whale{wn}.wloc(Iuse(i), :);

                    %calculate expected TDOA
                    R = sqrt(sum(hloc-wloc, 2).^2);

                    E.expTDOA_doa(ndet, 1) = (R(1)-R(2))./c;
                    E.expTDOA_doa(ndet, 2) = (R(1)-R(3))./c;
                    E.expTDOA_doa(ndet, 3) = (R(1)-R(4))./c;
                    E.expTDOA_doa(ndet, 4) = (R(2)-R(3))./c;
                    E.expTDOA_doa(ndet, 5) = (R(2)-R(4))./c;
                    E.expTDOA_doa(ndet, 6) = (R(3)-R(4))./c;

                    E.expTDOA_doa(ndet, 1) = (R(1)-R(2))./c;
                    E.expTDOA_doa(ndet, 2) = (R(1)-R(3))./c;
                    E.expTDOA_doa(ndet, 3) = (R(1)-R(4))./c;
                    E.expTDOA_doa(ndet, 4) = (R(2)-R(3))./c;
                    E.expTDOA_doa(ndet, 5) = (R(2)-R(4))./c;
                    E.expTDOA_doa(ndet, 6) = (R(3)-R(4))./c;

                    % calculate expected TDOA from model location:
                    R = sqrt(sum(hloc-whale{wn}.wloc(Iuse(i), :), 2).^2);
                    E.expTDOA_mod(ndet, 1) = (R(1)-R(2))./c;
                    E.expTDOA_mod(ndet, 2) = (R(1)-R(3))./c;
                    E.expTDOA_mod(ndet, 3) = (R(1)-R(4))./c;
                    E.expTDOA_mod(ndet, 4) = (R(2)-R(3))./c;
                    E.expTDOA_mod(ndet, 5) = (R(2)-R(4))./c;
                    E.expTDOA_mod(ndet, 6) = (R(3)-R(4))./c;

                    % calculate error from expected TDOA
                    E.measTDOA(ndet, :) = whale{wn}.TDOA(Iuse(i), 13:18);
                    E.TDet(ndet) = whale{wn}.TDet(Iuse(i));

                    for itdoa = 1:6
                        drft = polyval(dp{itdoa}, E.TDet(ndet));
                        E.drift(ndet, itdoa) = drft;
                        E.measTDOA_driftCorrected(ndet, itdoa) = E.measTDOA(ndet, itdoa) - drft;
                    end

                    %                     E.cest(ndet, 1) = (R(1) - R(2))./E.measTDOA_driftCorrected(ndet, 1);
                    %                     E.cest(ndet, 2) = (R(1) - R(3))./E.measTDOA_driftCorrected(ndet, 2);
                    %                     E.cest(ndet, 3) = (R(1) - R(4))./E.measTDOA_driftCorrected(ndet, 3);
                    %                     E.cest(ndet, 4) = (R(2) - R(3))./E.measTDOA_driftCorrected(ndet, 4);
                    %                     E.cest(ndet, 5) = (R(2) - R(4))./E.measTDOA_driftCorrected(ndet, 5);
                    %                     E.cest(ndet, 6) = (R(3) - R(4))./E.measTDOA_driftCorrected(ndet, 6);


                    %                     E.cest_nodrift(ndet, 1) = (R(1) - R(2))./E.measTDOA(ndet, 1);
                    %                     E.cest_nodrift(ndet, 2) = (R(1) - R(3))./E.measTDOA(ndet, 2);
                    %                     E.cest_nodrift(ndet, 3) = (R(1) - R(4))./E.measTDOA(ndet, 3);
                    %                     E.cest_nodrift(ndet, 4) = (R(2) - R(3))./E.measTDOA(ndet, 4);
                    %                     E.cest_nodrift(ndet, 5) = (R(2) - R(4))./E.measTDOA(ndet, 5);
                    %                     E.cest_nodrift(ndet, 6) = (R(3) - R(4))./E.measTDOA(ndet, 6);

                    % localize with large ap only
                    [~, Ibest] = min(sum((M.TDOA(:,13:18)-E.measTDOA_driftCorrected(ndet, :)).^2, 2));
                    E.wloc_lrg(ndet, :) = M.wloc(Ibest, :);
                end
            end
            ncluster = ncluster + 1;
            E.meanErr(ncluster, :) = mean(E.expTDOA_mod(nstart:ndet, :) - E.measTDOA(nstart:ndet, :));
            E.meanT(ncluster) = mean(E.TDet(nstart:ndet));
        end

    end
end

%%
figure(10)
for sp = 1:6
    subplot(6,1,sp)
    %     plot(E.TDet+0.2*randn(size(E.TDet)), (E.expTDOA_doa(:,sp) - E.measTDOA_driftCorrected(:,sp)), '.')
    %     ylim([0,.3])
    histogram(E.expTDOA_mod(:,sp) - (E.measTDOA(:,sp) -E.drift(:,sp)), -1:.01:1)
    hold on
    %     hold on
    %     histogram(E.expTDOA_doa(:,sp) - E.measTDOA(:,sp))
    %     hold off
    %     plot(E.expTDOA_doa(:,sp), 'o')
    %     hold on
    %     plot(E.measTDOA(:,sp), 'x')
    %     hold off
end
%%
figure(11)
for sp =1:6
    subplot(6,1,sp)
    plot(E.expTDOA_mod(:,sp) - E.expTDOA_doa(:,sp), '.')
    ylim([-0.000001, 0.000001])
end

figure(12)

for sp = 1:3
    subplot(3,1,sp)
    plot(E.wloc_doa(:,sp), 'x')
    hold on
    plot(E.wloc_mod(:,sp), '.')
    hold off
end

figure(13)
for sp = 1:6
    subplot(6,1,sp)
    plot(E.meanT, E.meanErr(:,sp), '.')
    grid on
    coef = polyfit(E.meanT, E.meanErr(:,sp), 1);
    errfit = polyval(coef, E.meanT);
    hold on
    plot(E.meanT, errfit)
    plot(tdrift, -drift(sp, :))
    hold off

end
legend('TDOAexp-TDOAmeas', 'linear fit', 'measured drift')

%%
