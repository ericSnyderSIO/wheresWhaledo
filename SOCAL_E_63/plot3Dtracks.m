load('D:\SOCAL_E_63\tracking\interns2022\ericEdits_allTracks\track600_180611_110414\SOCAL_E_63_track600_180611_110414_ericMod_localized_cleaned.mat')

load('D:\Writing\wheresWhaledo\figures\tracks\hydLoc.mat')

global brushing
loadParams('D:\MATLAB_addons\gitHub\wheresWhaledo\brushing.params')

load('D:\Writing\wheresWhaledo\figures\tracks\hydLoc.mat')

G = load('D:\SOCAL_E_63\bathymetry\SOCAL_E_63_GMRT.mat');

y2k = datenum([2000, 0, 0, 0, 0, 0]);

[x,~] = latlon2xy_wgs84(h0(1).*ones(size(G.lon)), G.lon, h0(1), h0(2));
[~,y] = latlon2xy_wgs84(G.lat, h0(2).*ones(size(G.lat)), h0(1), h0(2));

% define area of plot
gridLim = [-2100, 2100; -2100, 2100];
Ix = find(x>=-2100 & x<=2100);
Iy = find(y>=-2100 & y<=2100);

% 
% col = bone(100);
% col = col(10:50, :);

col = bone(1000);
% equation for color index (want color to be scaled as an inverted
% expontential):
I = -exp((1:100)./20);

% scale to be between 1 and 1000:
I = (I-min(I))/(max(I)-min(I))*(950-1) + 1;
I = round(fliplr(I));

col = col(I, :);

xv = gridLim(1,1):gridLim(1,2);
yv = gridLim(2,1):gridLim(2,2);


figure(63)

[cntr, hc] = contour3(x(Ix), y(Iy), G.z(Iy,Ix) - 20, -1340:10:-800, 'edgeColor', [.8, .8, .8]);
% clabel(cntr, hc, [-1340:20:-1250, -1200:100:-1000], 'color', [.6, .6, .6])
hold on
sc = surf(x(Ix), y(Iy), G.z(Iy,Ix)-20);
colormap(col)
shading interp

axis([-2000, 2000, -2000, 2000])
pbaspect([1,1,.5])
grid on
set(gca,'GridColor',[0.3 0.3 0.3])
xlabel('East to West (m)')
ylabel('North to South (m)')
zlabel('Depth (m)')
plot3(h(1:2,1), h(1:2,2), h(1:2, 3)-abs(h0(3)) + 6, 's', 'markerFaceColor', [.9,.9,.9], 'Color', [.5,.5,.5])
plot3(h(3:4,1), h(3:4,2), h(3:4, 3)-abs(h0(3)) + 10, 'o', 'markerFaceColor', [.9,.9,.9], 'Color', [.5,.5,.5])

% axis([-2000, 2000, -2000, 2000, -1350, -700])

for wn = 1:numel(whale)
    plot3(whale{wn}.wlocSmooth(:,1), whale{wn}.wlocSmooth(:,2), whale{wn}.wlocSmooth(:,3)-abs(h0(3)), ....
        'o', 'markerFacecolor', brushing.params.colorMat(wn+2, :), 'markerEdgeColor', 'none')
    
    % plot onto backplane
    scatter3(-2000.*ones(size(whale{wn}.wlocSmooth(:,1))), whale{wn}.wlocSmooth(:,2), whale{wn}.wlocSmooth(:,3)-abs(h0(3)), ...
        9, brushing.params.colorMat((wn+2).*ones(size(whale{wn}.TDet)), :), 'filled', 'MarkerFaceAlpha', 0.5)
end
hold off
