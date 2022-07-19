close all
clear all

map = imread("PIA22770.tif");
global size_map
size_map = size(map);
r0 = 2574.73e3;

%% plot the map
figure = tiledlayout(1,1, "Padding", "compact");

imshow(map)
equator_loc = lat_to_loc(0);
greenwich_loc = lon_to_loc(0);

%% plot final position of Huygens
hold on
%plot(lon_to_loc(185.43), lat_to_loc(-8.61), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
plot([lon_to_loc(185.43), lon_to_loc(196.0701)], [lat_to_loc(-8.61), lat_to_loc(-10.3213)], 'r--', 'LineWidth', 2)
plot(lon_to_loc(196.0701), lat_to_loc(-10.3213), 'r*', 'MarkerSize', 12, 'LineWidth', 2);
% landing site
plot(lon_to_loc(192.3247), lat_to_loc(-10.2507), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
hold off

%% integration until parachute
entryYelle = table2array(readtable("../output/entryYelle.csv"));
azim = 259.895; % positive North to East 
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
plot(lon_to_loc(174.2), lat_to_loc(-6.6), 'g^', 'MarkerSize', 12, 'LineWidth', 2);
hold off

hold on
plot([lon_to_loc(lon0-800e3*180/(pi*r0)), lon_to_loc(lon0+800e3*180/(pi*r0))], [lat_to_loc(lat0+800e3*180/(pi*r0)), lat_to_loc(lat0+800e3*180/(pi*r0))], 'g-.', 'MarkerSize', 12, 'LineWidth', 2)
plot([lon_to_loc(lon0-800e3*180/(pi*r0)), lon_to_loc(lon0+800e3*180/(pi*r0))], [lat_to_loc(lat0-800e3*180/(pi*r0)), lat_to_loc(lat0-800e3*180/(pi*r0))], 'g-.', 'MarkerSize', 12, 'LineWidth', 2)
plot([lon_to_loc(lon0-800e3*180/(pi*r0)), lon_to_loc(lon0-800e3*180/(pi*r0))], [lat_to_loc(lat0-800e3*180/(pi*r0)), lat_to_loc(lat0+800e3*180/(pi*r0))], 'g-.', 'MarkerSize', 12, 'LineWidth', 2)
plot([lon_to_loc(lon0+800e3*180/(pi*r0)), lon_to_loc(lon0+800e3*180/(pi*r0))], [lat_to_loc(lat0-800e3*180/(pi*r0)), lat_to_loc(lat0+800e3*180/(pi*r0))], 'g-.', 'MarkerSize', 12, 'LineWidth', 2)
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
       'reconstruction before parachute (integration)', 'parachute deployment (integration)', 'landing site (integration)', 'Selk', 'Antilia Faculae', 'Mindanao Facula', 'Field Of View' , FontSize=20);

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