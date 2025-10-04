clc
clear all
T0 = 518.7; % Sea-level standard temperature (R)
P0 = 2116.22; % Sea-level standard pressure (lb/ft^2)
rho0 = 0.0023769; % Sea-level standard density (slugs/ft^3)
g = 32.174;  % Gravity (ft/s^2)
R = 1716;    % Specific gas constant for air (ft^2/s^2 R)
gamma = 1.4; % Ratio of specific heats for air
C = 2.27e-7; % Sutherland constant (lb-s/ft^2 R^0.5)
S = 198.72;  % Sutherland's temperature constant (R)
a = -3.56616; % Lapse rate in troposphere (R/1000 ft)

% Altitudes (in feet)
altitudes = 0:1000:65000; % Altitudes from 0 to 65,000 ft in steps of 1000 ft
n = length(altitudes);

% Initialize result
T = zeros(1, n); % Temperature (R)
P = zeros(1, n); % Pressure (lb/ft^2)
rho = zeros(1, n); % Density (slugs/ft^3)
mu = zeros(1, n); % Dynamic viscosity (lb-s/ft^2)
T_ratio = zeros(1, n); % Temperature ratio (T/T0)
P_ratio = zeros(1, n); % Pressure ratio (P/P0)
rho_ratio = zeros(1, n); % Density ratio (rho/rho0)
P_psi = zeros(1, n); % Pressure in psi
Va = zeros(1, n); % Speed of sound (ft/s)

% Standard Atmosphere
for i = 1:n
    h = altitudes(i); % Current altitude
    if h <= 36000
        % Troposphere
        T(i) = T0 + a * h / 1000;
        P(i) = P0 * (T(i) / T0)^5.2561;
    else
        % Stratosphere
        T(i) = T(37); % Temperature at 37,000 ft
        P(i) = P(37)*exp(-((h - 36089) / (20806.7)));
    end
    % Density Calculation
    rho(i) = P(i)/(R*T(i));
    % Dynamic viscosity using Sutherland's equation
    mu(i) = C * T(i)^(1.5) / (T(i) + S);
    % Ratios
    T_ratio(i) = T(i) / T0;
    P_ratio(i) = P(i) / P0;
    rho_ratio(i) = rho(i) / rho0;
    % Convert Pressure to psi
    P_psi(i) = P(i) / 144; % 1 psi = 144 lb/ft^2
    % Speed of sound
    Va(i) = sqrt(gamma * R * T(i)); % Speed of sound (ft/s)
end

dataTable = table(altitudes', T', T_ratio', P_psi', P_ratio', rho', rho_ratio', mu', Va', ...
    'VariableNames', {'Altitude_ft', 'Temperature_R', 'Temp_Ratio', 'Pressure_psi', ...
    'Pressure_Ratio', 'Density_slugs_ft3', 'Density_Ratio', ...
    'Viscosity_lb_s_ft2', 'Speed_of_Sound_ft_s'});

disp(dataTable);
figure('Color','w','Name','US Std Atmosphere (0–65 kft)');
subplot(2,2,1); plot(T,altitudes/1000,'LineWidth',1.2);
grid on;
xlabel('T (R)');
ylabel('Altitude (kft)'); 
title('Temperature')
subplot(2,2,2);
plot(P_psi,altitudes/1000,'LineWidth',1.2);
grid on 
xlabel('P (psi)');
ylabel('Altitude (kft)');
title('Pressure')
subplot(2,2,3);
plot(rho,altitudes/1000,'LineWidth',1.2);
grid on 
xlabel('\rho (slug/ft^3)');
ylabel('Altitude (kft)');
title('Density')
subplot(2,2,4);
plot(Va,altitudes/1000,'LineWidth',1.2);
grid on 
xlabel('a (ft/s)');
ylabel('Altitude (kft)');
title('Speed of Sound')
sgtitle('US Standard Atmosphere (0–65 kft)')
Vtrue_kt = 300;
Vtrue    = Vtrue_kt * 1.68781;     % ft/s
Mach_at_Vtrue = Vtrue ./ Va;       % Mach vs altitude
