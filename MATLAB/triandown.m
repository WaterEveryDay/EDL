% Seb Henry Jul 22
map = imread("PIA22770.tif");
global size_map r0
size_map = size(map);
r0 = 2574.73e3;

disp("###############################")
disp("########## DOWN CASE ##########")
disp("###############################")
%% integration until ground
entryYelle = table2array(readtable("../output/entryYelle_nonplanar.csv"));
% lon0 = entryYelle(1,6)*180/pi;
% lat0 = entryYelle(1,7)*180/pi;
lonf = entryYelle(end,6)*180/pi; % ground
latf = entryYelle(end,7)*180/pi;

%% TRIANGULATION AT 1250KM ALTITUDE
lon0 = 185.43; % deg
lat0 = -8.61; % deg
alt0 = 1247.68717e3; % m
azim = 259.895; % deg
fph = 65.62; % deg
sigma = 90*pi/180 /(1024) / 10; % rad

% transformation to Titan
rsat_T = lonlatalt2Titan(lon0, lat0, alt0);

% transformation to the RSW
R = rsat_T/norm(rsat_T);
W = cross(R, [sind(-azim); cosd(azim); 0])/norm(cross(R, [sind(-azim); cosd(azim); 0]));
S = cross(W, R);

A_T_to_RSW = [R'; S'; W'];
A_RSW_to_C = [cosd(30), 0, sind(30); 0, 1, 0; -sind(30), 0, cosd(30)];
A_T_to_C = A_RSW_to_C * A_T_to_RSW;

% list of features
r_features_Titan = [
    lonlatalt2Titan(199, 7, 0)'; % Selk NOT VISIBLE with canted Camera
    lonlatalt2Titan(189, -11, 0)'; % Antilia
    lonlatalt2Titan(174.2, -6.6, 0)'; % Mindanao
    lonlatalt2Titan(164.1, -10.4, 0)'; % Shikoku
    lonlatalt2Titan(162, -27, 0)'; % Perkuna    
]';

%% generate covariances
Rxx = sigma^2*eye(2);

r_features_RSW = A_T_to_RSW * r_features_Titan;
rsat_RSW = A_T_to_RSW * rsat_T;

l_RSW = r_features_RSW - rsat_RSW;
l_C = A_RSW_to_C * l_RSW;
x = l_C ./ l_C(3,:);

A_list = cat(3,A_RSW_to_C, A_RSW_to_C);
Rxx_list = cat(3, Rxx, Rxx);
for i = 3:length(l_C(1,:))
    A_list = cat(3,A_list, A_RSW_to_C);
    Rxx_list = cat(3, Rxx_list, Rxx);
end

[r_hat_RSW_LOST, cov_LOST_init] = TriLOST2(x, r_features_RSW, A_list, Rxx_list);
[r_hat_RSW_DLT, cov_DLT_init] = TriDLT2(x, r_features_RSW, A_list, Rxx_list);

%% Prints
disp("### LOST ###")
disp("cov_LOST_init = ")
disp(cov_LOST_init)
disp("cov_LOST_end (h, lon, lat) = ")
cov_LOST_end_lonlat = [0.000922109 8.95829e-07 7.94637e-08 % from c++
                        8.95829e-07 1.49653e-09 5.87896e-11
                        7.94637e-08 5.87896e-11 6.22308e-10];
disp(cov_LOST_end_lonlat)
disp("sigma_r LOST init")
disp(sqrt(trace(cov_LOST_init)))

disp("ground ellipse features LOST")
[xlost, ylost] = getErrorEllipse([lonf, latf], cov_LOST_end_lonlat(2:3,2:3)*(180/pi)^2, 0.99);

disp("### DLT ###")
disp("cov_DLT_init = ")
disp(cov_DLT_init)

disp("cov_DLT_end (h, lon, lat) = ") % crom c++
cov_DLT_end_lonlat = [0.000925145 8.98546e-07 8.54181e-08
                      8.98546e-07 1.50598e-09 6.96766e-11
                      8.54181e-08 6.96766e-11 6.41662e-10];
disp(cov_DLT_end_lonlat)

disp("sigma_r DLT init")
disp(sqrt(trace(cov_DLT_init)))

disp("ground ellipse features DLT")
[xdlt, ydlt] = getErrorEllipse([lonf, latf], cov_DLT_end_lonlat(2:3,2:3)*(180/pi)^2, 0.99);

%% plot the map
%figure = tiledlayout(1,1, "Padding", "compact");
imshow(map)

hold on
plot(lon_to_loc(xlost), lat_to_loc(ylost), '-', 'Color','g', 'LineWidth',3)
plot(lon_to_loc(xdlt), lat_to_loc(ydlt), '--', 'Color','#5F8575', 'LineWidth',3)

mc_lost = table2array(readtable("../output/MC_LOST_entry_down.csv"));
x = lon_to_loc(mc_lost(:,2)*180/pi);
y = lat_to_loc(mc_lost(:,3)*180/pi);

scatter(x, y, '+', MarkerEdgeColor='g')

mc_dlt = table2array(readtable("../output/MC_DLT_entry_down.csv"));
x = lon_to_loc(mc_dlt(:,2)*180/pi);
y = lat_to_loc(mc_dlt(:,3)*180/pi);
scatter(x, y, '*', MarkerEdgeColor='#5F8575')

% landing site
plot(lon_to_loc(entryYelle(:,6)*180/pi), lat_to_loc(entryYelle(:,7)*180/pi), 'b-', 'LineWidth', 2)
plot(lon_to_loc(lonf), lat_to_loc(latf), 'b+', 'MarkerSize', 12, 'LineWidth', 2);

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
xticklabels(string(lon_ticks))
xlabel('longitude (°)', 'FontSize',20)

lat_ticks = 90:-0.01:-90;
yticks(lat_to_loc(lat_ticks))
yticklabels(string(lat_ticks))
ylabel('latitude (°)', 'FontSize',20)

legend('ellipse LOST', 'ellipse DLT', 'MC LOST', 'MC DLT', 'trajectory (integration)', 'landing site (integration)' , FontSize=20);

xlim(lon_to_loc([196.345, 196.31]))
ylim(lat_to_loc([-10.52, -10.54]))

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