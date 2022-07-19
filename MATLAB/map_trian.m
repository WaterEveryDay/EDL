close all
clear all

map = imread("PIA22770.tif");
global size_map r0
size_map = size(map);
r0 = 2574.73e3;

%% plot the map
figure = tiledlayout(1,1, "Padding", "compact");

imshow(map)

%%
lon0 = 185.43;
lat0 = -8.61;
alt0 = 1247.68717e3;
rsat_T = lonlatalt2Titan(lon0, lat0, alt0);
azim = 259.895;
R = rsat_T/norm(rsat_T);
W = cross(R, [sind(-azim); cosd(azim); 0])/norm(cross(R, [sind(-azim); cosd(azim); 0]));
S = cross(W, R);

A_T_to_RSW = [R'; S'; W'];
A_RSW_to_C = [0, 0, 1; 0, 1, 0; -1, 0, 0];

A_T_to_C = A_RSW_to_C * A_T_to_RSW;


rSelk_T = lonlatalt2Titan(199, 7, 0);
rAntilia_T = lonlatalt2Titan(189, -11, 0);
rMindanao_T = lonlatalt2Titan(174.2, 6.6, 0);

rsat_RSW = A_T_to_RSW * rsat_T;
rSelk_RSW = A_T_to_RSW * rSelk_T;
rAntilia_RSW = A_T_to_RSW *rAntilia_T;
rMindanao_RSW = A_T_to_RSW *rMindanao_T;

k = [0;0;1];
S = [1 0 0; 0 1 0];
l1_RSW = (rSelk_RSW - rsat_RSW);
l2_RSW = (rAntilia_RSW - rsat_RSW);
l3_RSW = (rMindanao_RSW - rsat_RSW);
l1_C = A_RSW_to_C * l1_RSW;
l2_C = A_RSW_to_C * l2_RSW;
l3_C = A_RSW_to_C * l3_RSW;
x1 = l1_C / (k'*l1_C);
x2 = l2_C / (k'*l2_C);
x3 = l3_C / (k'*l3_C);

sigma = 90*pi/180 / (1024) / 10;
R_xx = sigma^2*eye(2);
[r_hat_RSW, cov_LOST] = TriLOST2([x1, x2, x3], [rSelk_RSW, rAntilia_RSW, rMindanao_RSW], cat(3,A_RSW_to_C, A_RSW_to_C, A_RSW_to_C), cat(3, R_xx, R_xx, R_xx))
[r_hat_RSW, cov_DLT] = TriDLT2([x1, x2, x3], [rSelk_RSW, rAntilia_RSW, rMindanao_RSW], cat(3,A_RSW_to_C, A_RSW_to_C, A_RSW_to_C), cat(3, R_xx, R_xx, R_xx))
sqrt(trace(cov_LOST))
sqrt(trace(cov_DLT))

%% integration until parachute
entryYelle = table2array(readtable("../output/entryYelle.csv"));
 % positive North to East 
fph = 65.62;
lon0 = 185.43;
lat0 = -8.61;

dt = entryYelle(2,1) - entryYelle(1,1);
lonf = lon0;
latf = lat0;

for i = 1:2686
    h = entryYelle(i,4);
    dS = entryYelle(i,2).*cos(entryYelle(i,3)) * dt;
    dS = dS / (r0 + h) * r0;
    dx = dS*sind(-azim);
    dy = dS*cosd(azim);
    dlon = dx*360/(2*pi*r0) / cosd(latf);
    dlat = dy*180/(pi*r0);
    lonf = lonf+dlon;
    latf = latf+dlat;
end

hold on
%plot(lon_to_loc(lon0), lat_to_loc(lat0), 'b+', 'MarkerSize', 8, 'LineWidth', 2);
plot([lon_to_loc(lon0), lon_to_loc(lonf)], [lat_to_loc(lat0), lat_to_loc(latf)], 'b-', 'LineWidth', 2)
plot(lon_to_loc(lonf), lat_to_loc(latf), 'b*', 'MarkerSize', 12, 'LineWidth', 2);
hold off

lonf = lon0;
latf = lat0;

for i = 1:length(entryYelle(:,1))
    h = entryYelle(i,4);
    dS = entryYelle(i,2).*cos(entryYelle(i,3)) * dt;
    dS = dS / (r0 + h) * r0;
    dx = dS*sind(-azim);
    dy = dS*cosd(azim);
    dlon = dx*360/(2*pi*r0) / cosd(latf);
    dlat = dy*180/(pi*r0);
    lonf = lonf+dlon;
    latf = latf+dlat;
end
hold on
% landing site
plot(lon_to_loc(lonf), lat_to_loc(latf), 'b+', 'MarkerSize', 12, 'LineWidth', 2);
hold off

hold on
plot(lon_to_loc(199), lat_to_loc(7), 'go', 'MarkerSize', 12, 'LineWidth', 2);
plot(lon_to_loc(187), lat_to_loc(-11), 'gd', 'MarkerSize', 12, 'LineWidth', 2);
hold off


%% make the map pretty
axis on
ax = gca;
ax.FontSize=15;

grid on
ax = gca;
ax.GridColor = 'w';
ax.GridAlpha=1;

lon_ticks = 360:-5:0;
xticks(lon_to_loc(lon_ticks))
xticklabels(string(lon_ticks))
xlabel('longitude (°)', 'FontSize',20)

lat_ticks = 90:-5:-90;
yticks(lat_to_loc(lat_ticks))
yticklabels(string(lat_ticks))
ylabel('latitude (°)', 'FontSize',20)

legend('reconstruction before parachute (truth)', 'parachute deployment (truth)', 'landing site (truth)', ...
       'reconstruction before parachute (integration)', 'parachute deployment (integration)', 'landing site (integration)', 'Selk', 'Antilia Faculae' , FontSize=20);

scale = 0.02;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)


%% functions
function loc = lon_to_loc(lon_in_deg)
    global size_map
    size_map_2 = size_map(2);
    a = (1-size_map_2)/360;
    b = size_map_2;
    loc = a * lon_in_deg + b;
end

function loc = lat_to_loc(lat_in_deg)
    global size_map
    size_map_1 = size_map(1);
    a = (1-size_map_1)/180;
    b = size_map_1 /2 + 1/2;
    loc = a * lat_in_deg + b;
end

function r = lonlatalt2Titan(lon, lat, alt)
% angles in degrees
    global r0
    r = (r0+alt) * [cosd(lon)*cosd(lat); -sind(lon)*cosd(lat); sind(lat)];
end