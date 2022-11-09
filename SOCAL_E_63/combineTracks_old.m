global CTCparam
loadParams('CTCparams.txt')

% xwav tables:
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EE_C4_xwavLookupTable');
XH{1} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EW_C4_xwavLookupTable');
XH{2} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EN_xwavLookupTable');
XH{3} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_ES_xwavLookupTable');
XH{4} = xwavTable;

spd = 60*60*24; % seconds per day (converting from datenum)
numInst = 4; % number of instruments

df = dir('D:\SOCAL_E_63\tracking\interns2022\ericEdits_allTracks\*track*'); % directory of folders containing files

for ndf = 1:numel(df)

    d = dir(fullfile(df(ndf).folder, [df(ndf).name, '\*.mat']));
    if ~isempty(d)
        load(fullfile(d.folder, d.name))

        tstart = min([DET{1}.TDet(1), DET{2}.TDet(1)]);
        tend = max([DET{1}.TDet(end), DET{2}.TDet(end)]);

        % rerun detector with better lockout time:
        %         DET{3} = detectClicks_1ch(tstart, tend, XH{3}, 'D:\MATLAB_addons\gitHub\wheresWhaledo\detClicks_1ch.params');
        %         DET{4} = detectClicks_1ch(tstart, tend, XH{4}, 'D:\MATLAB_addons\gitHub\wheresWhaledo\detClicks_1ch.params');

        % click train corr to associate clicks on single channels

        whaleNums1 = unique(DET{1}.color);
        whaleNums2 = unique(DET{2}.color);

        whaleNums = unique([whaleNums1; whaleNums1]);

        for wn = 1:length(whaleNums)

            I1 = find(DET{1}.color==whaleNums(wn));
            whale{whaleNums(wn)-2} = DET{1}(I1, :);
            whale{1}.ctcTDOA = nan(length(I1), 6);
            for i1 = 1:length(I1)
                % define bounds of click train:
                tstart = DET{1}.TDet(I1(i1)) - CTCparam.twin/(2*spd);
                tend = tstart + CTCparam.twin/spd;

                % define time vector:
                tc = tstart:1/(spd*CTCparam.fsct):tend;

                X = zeros(length(tc), numInst);

                for ninst = 1:2
                    Idet = find(DET{ninst}.TDet>=tstart & DET{ninst}.TDet<=tend & DET{ninst}.color==whaleNums(wn));
                    for idet = 1:length(Idet)
                        % find index in tc which most closely matches
                        % detection time:
                        [~, I] = min((tc-DET{ninst}.TDet(Idet(idet))).^2);
                        X(I, ninst) = 1;
                    end
                    Xhann(:, ninst) = conv(X(:, ninst), CTCparam.Wk);
                end

                for ninst = 3:numInst
                    Idet = find(DET{ninst}.TDet>=tstart & DET{ninst}.TDet<=tend);
                    for idet = 1:length(Idet)
                        % find index in tc which most closely matches
                        % detection time:
                        [~, I] = min((tc-DET{ninst}.TDet(Idet(idet))).^2);
                        X(I, ninst) = 1;
                    end
                end
                Xhann(:, ninst) = conv(X(:, ninst), CTCparam.Wk);

                % need a method that combines whaleAss with CTC, with
                % outputs:
                % DET, where DET{3} and DET{4} now have detections
                % labeled as correct whales according the whaleAss
                % (from both DET{1} labels and DET{2} labels)
                % whale{wn}.TDOA which has taken the clicks found on
                % each instrument, determined whether they are
                % associated, and calculated an approx. TDOA from the
                % CTC.

                % steps:
                % 1. xcorr click trains
                % 2. Associate detections, as done in whaleAss
                % 3. If a click from another inst is aligned with the
                % center click from current array:
                % 3a. save whale{wn}.TDOA
                % 3b. label DET{inst}.color = wn+2

                tdoa = nan(3,1);
                for ninst = 2:4 % iterate through other instruments
                    [xctc, lags] = xcorr(Xhann(:, 1), Xhann(:, ninst), CTCparam.maxLag); % xcorr click trains

                    [pks, locs] = findpeaks(xctc, 'NPeaks', 2, 'SortStr', 'descend');
                    if length(pks)==2
                        if pks(1)*CTCparam.peakRatio>pks(2) % click trains correlated enough to assume they are the same whale

                            % TDOA estimate from CTC:
                            tdoa(ninst-1) = lags(locs(1))/CTCparam.fsct;

                            % Associate click
                            [terr, Ibest] = min(abs(DET{ninst}.TDet - DET{1}.TDet(I1(i1)) + tdoa(ninst-1)/spd));

                            if terr<=(CTCparam.tCloseEnough/spd)
                                DET{ninst}.Label(Ibest) = num2str(whaleNums(wn)-2);
                                DEt{ninst}.color(Ibest) = whaleNums(wn);
                            end

                        end
                    end
                end

                whale{whaleNums(wn)-2}.ctcTDOA_a1(i1, 1:3) = tdoa;
%                 whale{whaleNums(wn)-2}.ctcTDOA(i1, 4) = tdoa(2)-tdoa(1);
%                 whale{whaleNums(wn)-2}.ctcTDOA(i1, 5) = tdoa(3)-tdoa(1);
%                 whale{whaleNums(wn)-2}.ctcTDOA(i1, 6) = tdoa(3)-tdoa(2);

            end
            ok = 1;

        end
    end



end


