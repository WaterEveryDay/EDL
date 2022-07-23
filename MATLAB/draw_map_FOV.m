close all
clear all

map = imread("PIA22770.tif");
global size_map
size_map = size(map);
r0 = 2574.73e3;

%% plot the map
figure = tiledlayout(1,1, "Padding", "compact");

imshow(map)

%% integration until parachute
entryYelle = table2array(readtable("../output/entryYelle_nonplanar.csv"));
azim = 259.895; % positive North to East 
fph = 65.62;
lon0 = entryYelle(1,6)*180/pi;
lat0 = entryYelle(1,7)*180/pi;
lonp = entryYelle(2686,6)*180/pi; % parachute deployment
latp = entryYelle(2686,7)*180/pi;
lonf = entryYelle(end,6)*180/pi; % ground
latf = entryYelle(end,7)*180/pi;

lon800 = entryYelle(823,6)*180/pi;
lat800 = entryYelle(823,7)*180/pi;

hold on
plot(lon_to_loc(entryYelle(:,6)*180/pi), lat_to_loc(entryYelle(:,7)*180/pi), 'b-', 'LineWidth', 2)
%plot(lon_to_loc(entryYelle(1:2686,6)*180/pi), lat_to_loc(entryYelle(1:2686,7)*180/pi), 'b-', 'LineWidth', 2)
% landing site
plot(lon_to_loc(lonf), lat_to_loc(latf), 'b+', 'MarkerSize', 12, 'LineWidth', 2);
hold off

hold on
plot(lon_to_loc(199), lat_to_loc(7), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
plot(lon_to_loc(187), lat_to_loc(-11), 'rs', 'MarkerSize', 12, 'LineWidth', 2);
plot(lon_to_loc(174.2), lat_to_loc(-6.6), 'rd', 'MarkerSize', 12, 'LineWidth', 2);
plot(lon_to_loc(164.1), lat_to_loc(-10.4), 'r^', 'MarkerSize', 12, 'LineWidth', 2);
plot(lon_to_loc(162), lat_to_loc(-27), 'r>', 'MarkerSize', 12, 'LineWidth', 2);
hold off

hold on
%% green FOV
FOV1 = [-1200e3*180/(pi*r0), +1200e3*180/(pi*r0), +1200e3*180/(pi*r0), -1200e3*180/(pi*r0), -1200e3*180/(pi*r0);
    -1200e3*180/(pi*r0), -1200e3*180/(pi*r0), +1200e3*180/(pi*r0), +1200e3*180/(pi*r0), -1200e3*180/(pi*r0)];
rot_FOV = [cosd(azim+90), -sind(azim+90); sind(azim+90), cosd(azim+90)];
FOV1 = rot_FOV * (FOV1);
FOV1 = FOV1 + [lon0; lat0];
plot(lon_to_loc(FOV1(1,:)), lat_to_loc(FOV1(2,:)), 'g-', 'MarkerSize', 12, 'LineWidth', 2)

%% ORANGE FOV
% A_RSW_to_C = [cosd(45), 0, sind(45); 0, 1, 0; -sind(45), 0, cosd(45)];
FOV2 = [-1200e3*180/(pi*r0), +1200e3*180/(pi*r0), +1200e3*180/(pi*r0), -1200e3*180/(pi*r0), -1200e3*180/(pi*r0);
    -1200e3*tand(60)*180/(pi*r0), -1200e3*tand(60)*180/(pi*r0), 1200e3*tand(30)*180/(pi*r0), 1200e3*tand(30)*180/(pi*r0), -1200e3*tand(60)*180/(pi*r0)];
rot_FOV = [cosd(azim+90), -sind(azim+90); sind(azim+90), cosd(azim+90)];
FOV2 = rot_FOV * (FOV2);
FOV2 = FOV2 + [lon0; lat0];
plot(lon_to_loc(FOV2(1,:)), lat_to_loc(FOV2(2,:)), '--', 'Color','#ffa500', 'MarkerSize', 12, 'LineWidth', 2)
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

legend('reconstruction (integration)', 'landing site (integration)', ...
       'Selk Crater', 'Antilia Faculae', 'Mindanao Facula', 'Shikoku Facula', 'Perkunas Virgae',  ...
       'Field Of View straight down', 'Field Of View at 30° canted', FontSize=20);

scale = 0.02;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)
xlim(lon_to_loc([220, 120]))
ylim(lat_to_loc([25, -60]))


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