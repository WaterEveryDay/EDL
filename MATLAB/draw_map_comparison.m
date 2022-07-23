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
plot(lon_to_loc(entryYelle(1:2686,6)*180/pi), lat_to_loc(entryYelle(1:2686,7)*180/pi), 'b-', 'LineWidth', 2)
plot(lon_to_loc(lonp), lat_to_loc(latp), 'b*', 'MarkerSize', 12, 'LineWidth', 2);
% landing site
plot(lon_to_loc(lonf), lat_to_loc(latf), 'b+', 'MarkerSize', 12, 'LineWidth', 2);
hold off

%% plot final position of Huygens
entryTrue = table2array(readtable("../output/reconstructed_data.csv", 'VariableNamingRule','preserve'));
hold on
plot(lon_to_loc(entryTrue(:,4)), lat_to_loc(entryTrue(:,5)), 'r--', 'LineWidth', 2);
plot(lon_to_loc(196.0701), lat_to_loc(-10.3213), 'r*', 'MarkerSize', 12, 'LineWidth', 2);
% landing site
plot(lon_to_loc(192.3247), lat_to_loc(-10.2507), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
hold off

%% make the map pretty
axis on
ax = gca;
ax.FontSize=15;

grid on
ax = gca;
ax.GridColor = 'w';
ax.GridAlpha=1;

lon_ticks = 360:-2:0;
xticks(lon_to_loc(lon_ticks))
xticklabels(string(lon_ticks))
xlabel('longitude (°)', 'FontSize',20)

lat_ticks = 90:-2:-90;
yticks(lat_to_loc(lat_ticks))
yticklabels(string(lat_ticks))
ylabel('latitude (°)', 'FontSize',20)

legend('reconstruction before parachute (integration)', 'parachute deployment (integration)', 'landing site (integration)', ...
       'reconstruction before parachute (Kazeminejad et al.)', 'parachute deployment (Kazeminejad et al.)', 'landing site (Kazeminejad et al.)', ...
       FontSize=20, location='best');

scale = 0.01;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
pos(1) = pos(1)+scale*pos(3);
pos(3) = (1-scale)*pos(3);
set(gca, 'Position', pos)
xlim(lon_to_loc([198, 185]))
ylim(lat_to_loc([-7, -12]))


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