% Seb Henry Jul 22
close all 
clear all

global size_map r0
r0 = 2574.73e3;
map = imread("PIA22770.tif");
size_map = size(map);

%% integration until ground
entryYelle = table2array(readtable("../output/entryYelle_nonplanar.csv"));
lonf = entryYelle(end,6)*180/pi; % ground
latf = entryYelle(end,7)*180/pi;
imshow(map)
hold on
plot(lon_to_loc(entryYelle(:,6)*180/pi), lat_to_loc(entryYelle(:,7)*180/pi), 'b-', 'LineWidth', 2)
hold off


%% 
C_R = [3.21469e-10 -5.67693e-11; -5.67693e-11  1.00251e-11];
C_S = [6.85705e-10, -1.20569e-10 ; -1.20569e-10, 2.12e-11];
C_W = [1.50154e-11, 1.00112e-10; 1.00112e-10, 6.67476e-10];
hold on
disp("### R ###")
[x, y] = getErrorEllipse([lonf, latf], C_R*(180/pi)^2, 0.99);
x = lon_to_loc(x);
y = lat_to_loc(y);
disp("### S ###")
plot(x,y, LineWidth=2.5, Color='r')
[x, y] = getErrorEllipse([lonf, latf], C_S*(180/pi)^2, 0.99);
x = lon_to_loc(x);
y = lat_to_loc(y);
plot(x,y, '.', LineWidth=2, Color='g')
disp("### W ###")
[x, y] = getErrorEllipse([lonf, latf], C_W*(180/pi)^2, 0.99);
x = lon_to_loc(x);
y = lat_to_loc(y);
plot(x,y, '-.', LineWidth=2, Color='#ffa500')
hold off


%% make the map pretty
axis on
ax = gca;
ax.FontSize=15;

grid on
ax = gca;
ax.GridColor = 'w';
ax.GridAlpha=1;

lon_ticks = 360:-0.01:0;
xticks(lon_to_loc(lon_ticks))
xticklabels(string(lon_ticks*pi/180*r0))
xlabel('longitude (m)', 'FontSize',20)

lat_ticks = 90:-0.01:-90;
yticks(lat_to_loc(lat_ticks))
yticklabels(string(lat_ticks*pi/180*r0))
ylabel('latitude (m)', 'FontSize',20)

legend("reconstructed trajectory (integration)",  "ellipsoid pure altitude", ...
    "ellipsoid pure along-track", "ellipsoid pure cross-track", FontSize=20);

scale = 0.02;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)

xlim(lon_to_loc([196.335, 196.31]))
ylim(lat_to_loc([-10.52, -10.54]))

function [x, y] = getErrorEllipse(mu, Sigma, p)
% https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
    global r0;
    s = -2 * log(1 - p);
    [V, D] = eig(Sigma * s);

    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];

    sMa_angle = sqrt(max(abs(diag(D)))) ;
    sma_angle = sqrt(min(abs(diag(D)))) ;
    sMa = sMa_angle*pi/180*r0;
    disp("semi major axis = " + sMa + 'm')
    sma = sma_angle*pi/180*r0;
    disp("semi minor axis = " + sma + 'm')
    disp("standard deviation = " + sqrt(sum(diag(D))/s)*pi/180*r0 + "m")

    x = a(1, :) + mu(1);
    y = a(2, :) + mu(2);
end

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
