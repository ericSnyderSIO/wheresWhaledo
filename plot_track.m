function fig = plot_track()
% takes labeled points in figure 139, and estimates the whale's location
% based off of these points. Plots these estimates on a 3D plot, and
global REMORA

col = [0.1, 0.4, 0.8;
    0.9, 0.1, 0.3;
    0.25, 0.5, 0.05;
    0.7, 0.1, 0.6;
    0.9, 0.4, 0.1;
    0.27, 0.94, 0.94;
    0.1, 0.9, 0.1;
    0.9, 0.74, 0.95;
    0, 0, 0];

ms = 6;

spd = 24*60*60; % seconds per day (to convert between datenum and seconds

H1 = REMORA.H1;
H2 = REMORA.H2;
H0 = REMORA.H0;

H1(3) = -abs(H1(3)); % make depth negative
H2(3) = -abs(H2(3));
H0(3) = -abs(H0(3));

[h1(1), h1(2)] = latlon2xy(H1(1), H1(2), H0(1), H0(2)); % set hyd 1 loc in meters
h1(3) = H0(3) - H1(3); % set h1 depth
[h2(1), h2(2)] = latlon2xy(H2(1), H2(2), H0(1), H0(2)); % set hyd 2 loc in meters
h2(3) = H0(3) - H2(3);

% lsqnonlin options:
% options = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective');
% options.Display = 'none';

numwhale = 0;
for wn = 1:8 % iterate through each whale number
    
    % find points labeled as whale 'wn'
    N1 = find(REMORA.brushing.all.AR1.label==wn);
    N2 = find(REMORA.brushing.all.AR2.label==wn);
    
    REMORA.track{wn} = [];
    
    if ~isempty(N1)&&~isempty(N2)
        fprintf(['Calculating Whale Track: whale #', num2str(wn), '\n'])
        
        % Array 1:
        AZ1 = REMORA.brushing.all.AR1.Ang(N1,1);
        EL1 = REMORA.brushing.all.AR1.Ang(N1,2);
        T1 = REMORA.brushing.all.AR1.TDet(N1);
        
        % Array 2:
        AZ2 = REMORA.brushing.all.AR2.Ang(N2,1);
        EL2 = REMORA.brushing.all.AR2.Ang(N2,2);
        T2 = REMORA.brushing.all.AR2.TDet(N2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % smooth az/el
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % interpolate az/el
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % create interp time vector using only period where detections
        % overlap
        t = max([T1(1), T2(1)]):1/spd:min([T1(end), T2(end)]);
%         t = unique(sort([T1; T2])).';

        % interpolate. 'unwrap' included to correct for when whale crosses
        % 0/360 azimuth.
        az1 = interp1(T1, (unwrap(AZ1.*pi/180).*180/pi), t);
        el1 = interp1(T1, (EL1), t);
        az2 = interp1(T2, (unwrap(AZ2.*pi/180).*180/pi), t);
        el2 = interp1(T2, (EL2), t);
        
        REMORA.track{wn}.AR1.Ang = [AZ1, EL1];
        REMORA.track{wn}.AR1.t = T1;
        REMORA.track{wn}.AR2.Ang = [AZ2, EL2];
        REMORA.track{wn}.AR2.t = T2;
        
        REMORA.track{wn}.AR1.AngInterp = [az1; el1].';
        REMORA.track{wn}.AR1.tInterp = t.';
        
        REMORA.track{wn}.AR2.AngInterp = [az2; el2].';
        REMORA.track{wn}.AR2.tInterp = t.';
      
        N = length(t); % number of detections/samples/animal locations

        x = zeros(N,1);
        y = zeros(N,1);
        z = zeros(N,1);
        dL = zeros(N,1);
        skewDist = zeros(N,1);
        Lerr = zeros(N,3);
%         skewDist{wn} = zeros(N,1);
        
%         whaleLoc{wn} = zeros(N,3);
%         resnorm{wn} = zeros(N,1);
%         residual{wn} = zeros(N,3);
        %     for k = 1:10   % for testing
        for k = 1:N
            s1 = [sind(180-el1(k))*cosd(az1(k)), sind(180-el1(k))*sind(az1(k)), cosd(180-el1(k))];
            s2 = [sind(180-el2(k))*cosd(az2(k)), sind(180-el2(k))*sind(az2(k)), cosd(180-el2(k))];
            
            S = [s1; -s2].';
            dh = (h2-h1).';
            
            R = S\dh;
            
            L1 = h1+R(1)*s1;
            L2 = h2+R(2)*s2;
            
            whaleLoc = (L1+L2)./2;
           
            
            dL(k) = sqrt(sum((L1-L2).^2));
            skewDist(k) = sqrt(sum(((h1-h2)*(cross(s1, s2).')./sqrt(sum(cross(s1,s2).^2))).^2));
            Lerr(k,:) = L1-L2;
            
            x(k) = whaleLoc(1);
            y(k) = whaleLoc(2);
            z(k) = whaleLoc(3);
            
        end % for k
        
        [lat, lon] = xy2latlon(x, y, H0(1), H0(2));
        
        REMORA.track{wn}.whaleLoc = [x, y, abs(H0(3))-z];
        REMORA.track{wn}.whaleLatLon = [lat, lon, abs(H0(3))-z];
        REMORA.track{wn}.Lerr = Lerr;
        REMORA.track{wn}.skewDist = skewDist;
        REMORA.track{wn}.t = t;
 
        numwhale = numwhale+1;
    end
    
end



%% Generate plots

figno = 140;

% initialize axes values:
xmin = min([h1(1), h2(1)]);
xmax = max([h1(1), h2(1)]);
ymin = min([h1(2), h2(2)]);
ymax = max([h1(2), h2(2)]);
zmin = max([abs(H1(3)), abs(H2(3))]);
zmax = max([h1(3), h2(3)]);

fig = figure(figno);
pltwn = 0; % plotted whale number
for wn = 1:8
    if ~isempty(REMORA.track{wn})
        pltwn = pltwn+1;
        
        colscale{wn}(:,1) = linspace(max([col(wn,1).*0.5, 0]), min([col(wn,1).*1.5, 1]), length(REMORA.track{wn}.t));
        colscale{wn}(:,2) = linspace(max([col(wn,2).*0.5, 0]), min([col(wn,2).*1.5, 1]), length(REMORA.track{wn}.t));
        colscale{wn}(:,3) = linspace(max([col(wn,3).*0.5, 0]), min([col(wn,3).*1.5, 1]), length(REMORA.track{wn}.t));
        
        scatter3(REMORA.track{wn}.whaleLoc(:,1), REMORA.track{wn}.whaleLoc(:,2), REMORA.track{wn}.whaleLoc(:,3),  ms, colscale{wn}, 'filled')
        hold on
%         leg(pltwn) = {['Whale ', num2str(wn)]};
        
        % update maxima and minima for axes:
        xmin = min([REMORA.track{wn}.whaleLoc(:,1); xmin]);
        xmax = max([REMORA.track{wn}.whaleLoc(:,1); xmax]);
        ymin = min([REMORA.track{wn}.whaleLoc(:,2); ymin]);
        ymax = max([REMORA.track{wn}.whaleLoc(:,2); ymax]);
        zmin = min([REMORA.track{wn}.whaleLoc(:,3); zmin]);
        zmax = max([REMORA.track{wn}.whaleLoc(:,3); zmax]);
        
    end
end

scatter3(h1(1), h1(2), -H0(3), 'r^')
scatter3(h2(1), h2(2), -H0(3), 'r^')
hold off

% flip z axis
set(gca, 'zdir', 'reverse')

% set axes to be cubic:
maxaxis = max([xmax-xmin, ymax-ymin, zmax-zmin]);
xmean = (xmax+xmin)/2;
ymean = (ymax+ymin)/2;
zmean = (zmax+zmin)/2;

xlim([xmean - maxaxis/2, xmean + maxaxis/2])
ylim([ymean - maxaxis/2, ymean + maxaxis/2])
zlim([zmean - maxaxis/2, zmean + maxaxis/2])
pbaspect([1,1,1])

xlabel('E-W (m)')
ylabel('N-S (m)')
zlabel('Depth (m)')
title('3-D Whale Track Estimates')
% legend(leg)

figno = figno+1;
fig = figure(figno);
spn = 0; % subplot number
for wn = 1:8
    if ~isempty(REMORA.track{wn})
        if ~isempty(REMORA.track{wn}.t)
            spn = spn+1;
            subplot(numwhale,1,spn)
            
            
            fp = fill([REMORA.track{wn}.t, fliplr(REMORA.track{wn}.t)], ...
                [REMORA.track{wn}.whaleLoc(:,3) - abs(REMORA.track{wn}.skewDist./2); ...
                flipud(REMORA.track{wn}.whaleLoc(:,3) + abs(REMORA.track{wn}.skewDist./2))], [0.8, 0.8, 0.8]);
            fp.EdgeColor = 'none';
            hold on
            
            scatter(REMORA.track{wn}.t, REMORA.track{wn}.whaleLoc(:,3), ms, colscale{wn}, 'filled')
            hold off
            axis ij
            datetick
            grid on
            legend('Distance between DOA beams', 'Estimated whale depth')
            ylabel('Depth (m)')
            title(['Depth vs. Time, whale #', num2str(wn)])
        end
    end
end

xlabel('Time')

% figno = figno+1;
% fig = figure(figno);
% for wn = 1:8
%     if ~isempty(REMORA.track{wn})
%         if ~isempty(REMORA.track{wn}.t)
%             spn = spn+1;
%             clear colscale
%             colscale(:,1) = linspace(max([col(wn,1).*0.5, 0]), min([col(wn,1).*1.5, 1]), length(REMORA.track{wn}.t));
%             colscale(:,2) = linspace(max([col(wn,2).*0.5, 0]), min([col(wn,2).*1.5, 1]), length(REMORA.track{wn}.t));
%             colscale(:,3) = linspace(max([col(wn,3).*0.5, 0]), min([col(wn,3).*1.5, 1]), length(REMORA.track{wn}.t));
%             
%             scatter(REMORA.track{wn}.whaleLoc(:,1), REMORA.track{wn}.whaleLoc(:,2), ms, colscale(wn,:), 'filled')
%             hold on
%             plot(REMORA.track{wn}.whaleLoc(:,1)-abs(REMORA.track{wn}.Lerr(:,1)./2), ...
%                 REMORA.track{wn}.whaleLoc(:,2)-abs(REMORA.track{wn}.Lerr(:,2)./2), 'k:')
%             plot(REMORA.track{wn}.whaleLoc(:,1)+abs(REMORA.track{wn}.Lerr(:,1)./2), ...
%                 REMORA.track{wn}.whaleLoc(:,2)+abs(REMORA.track{wn}.Lerr(:,2)./2), 'k:')
%             
%         end
%     end
% end
% hold off
% title('x-y plot')
% xlabel('x (m)')
% ylabel('y (m)')
% legend('Estimated whale location', 'Distance between DOA beams')
% grid on

