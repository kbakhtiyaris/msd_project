% =========================================================================
% SMART PLANT WATERING SYSTEM - MATLAB SIMULATION
% Multi-Zone Irrigation with Adaptive Control
% Mechatronics System Design Project
% =========================================================================

clear all; close all; clc;

%% SIMULATION PARAMETERS
sim_time = 24*7;          % Simulation duration (hours) - 1 week
dt = 0.25;                % Time step (hours) - 15 minutes
time = 0:dt:sim_time;     % Time vector
num_zones = 4;            % Number of irrigation zones

%% PLANT PROFILES - Different plant types for each zone
% [min_moisture, max_moisture, optimal_temp, max_temp, min_humidity, max_humidity]
plant_profiles = {
    'Succulent',    [20, 40, 15, 30, 20, 40];   % Zone 1
    'Herbs',        [35, 55, 18, 26, 40, 60];   % Zone 2
    'Vegetables',   [45, 70, 18, 28, 50, 70];   % Zone 3
    'Ferns',        [50, 70, 18, 24, 60, 80];   % Zone 4
};

%% SYSTEM PARAMETERS
pump_flow_rate = 2.0;      % Liters per minute
irrigation_duration = 10;   % Minutes per watering event
zone_area = 0.5;           % Square meters per zone
soil_depth = 0.2;          % Meters
soil_volume = zone_area * soil_depth * 1000; % Liters

% Water added per irrigation event (liters)
water_per_event = (pump_flow_rate * irrigation_duration) / num_zones;

% Moisture increase per liter of water (%)
moisture_increase_per_liter = 100 / soil_volume;

%% ENVIRONMENTAL CONDITIONS (Simulated daily cycles)
% Temperature varies from 15°C (night) to 30°C (day)
temp_ambient = @(t) 22.5 + 7.5*sin(2*pi*t/24 - pi/2);

% Humidity varies from 40% (day) to 70% (night)
humidity_ambient = @(t) 55 - 15*sin(2*pi*t/24 - pi/2);

%% INITIALIZE VARIABLES
soil_moisture = zeros(num_zones, length(time));
temperature = zeros(1, length(time));
humidity = zeros(1, length(time));
pump_status = zeros(num_zones, length(time));
water_consumed = zeros(num_zones, length(time));

% Initial soil moisture (%)
soil_moisture(:,1) = [30; 40; 50; 55];

%% EVAPOTRANSPIRATION MODEL
% ET = Kc * ET0
% ET0 depends on temperature and humidity
function et_rate = calculate_ET(temp, humid, plant_type)
    % Base ET0 calculation (simplified Penman-Monteith)
    % Higher temperature and lower humidity increase ET
    et0 = 0.5 * (temp/25) * (1 - humid/100);
    
    % Crop coefficient based on plant type
    switch plant_type
        case 'Succulent'
            kc = 0.3;  % Low water loss
        case 'Herbs'
            kc = 0.7;  % Moderate water loss
        case 'Vegetables'
            kc = 0.9;  % High water loss
        case 'Ferns'
            kc = 1.0;  % Very high water loss
        otherwise
            kc = 0.7;
    end
    
    % ET rate as % moisture loss per hour
    et_rate = kc * et0;
end

%% DRAINAGE MODEL
function drain_rate = calculate_drainage(moisture)
    % Drainage increases when soil moisture exceeds field capacity (60%)
    if moisture > 60
        drain_rate = 0.5 * (moisture - 60) / 40;  % % per hour
    else
        drain_rate = 0;
    end
end

%% IRRIGATION DECISION LOGIC
function should_water = check_irrigation_needed(moisture, temp, humid, profile)
    % Unpack plant profile
    min_moist = profile(1);
    max_moist = profile(2);
    opt_temp = profile(3);
    max_temp = profile(4);
    min_hum = profile(5);
    max_hum = profile(6);
    
    % Basic threshold check
    below_threshold = moisture < min_moist;
    
    % Environmental stress factor
    temp_stress = (temp > max_temp);
    humid_stress = (humid < min_hum);
    
    % Adaptive threshold - water earlier if stressed
    if temp_stress || humid_stress
        adaptive_threshold = min_moist + 5;  % Water 5% earlier
    else
        adaptive_threshold = min_moist;
    end
    
    should_water = (moisture < adaptive_threshold);
end

%% MAIN SIMULATION LOOP
fprintf('Starting Smart Plant Watering System Simulation...\n');
fprintf('Duration: %d hours (%.1f days)\n\n', sim_time, sim_time/24);

for i = 2:length(time)
    % Update environmental conditions
    temperature(i) = temp_ambient(time(i));
    humidity(i) = humidity_ambient(time(i));
    
    % Process each zone
    for zone = 1:num_zones
        % Get plant profile
        plant_name = plant_profiles{zone, 1};
        profile = plant_profiles{zone, 2};
        
        % Calculate water losses
        et_loss = calculate_ET(temperature(i), humidity(i), plant_name);
        drain_loss = calculate_drainage(soil_moisture(zone, i-1));
        
        % Natural moisture decrease
        moisture_loss = et_loss + drain_loss;
        soil_moisture(zone, i) = soil_moisture(zone, i-1) - moisture_loss * dt;
        
        % Check if irrigation needed
        if check_irrigation_needed(soil_moisture(zone, i), temperature(i), ...
                                    humidity(i), profile)
            % Activate pump for this zone
            pump_status(zone, i) = 1;
            
            % Add water
            water_added = water_per_event;
            moisture_gain = water_added * moisture_increase_per_liter;
            soil_moisture(zone, i) = soil_moisture(zone, i) + moisture_gain;
            
            % Track water consumption
            water_consumed(zone, i) = water_added;
            
            % Ensure moisture doesn't exceed 100%
            if soil_moisture(zone, i) > 100
                soil_moisture(zone, i) = 100;
            end
        end
        
        % Ensure moisture doesn't go below 0%
        if soil_moisture(zone, i) < 0
            soil_moisture(zone, i) = 0;
        end
    end
end

fprintf('Simulation completed!\n\n');

%% CALCULATE STATISTICS
total_water_used = sum(water_consumed, 'all');
irrigation_events = sum(pump_status, 'all');

fprintf('=== SIMULATION RESULTS ===\n');
fprintf('Total water consumed: %.2f liters\n', total_water_used);
fprintf('Total irrigation events: %d\n', irrigation_events);
fprintf('Average water per zone: %.2f liters\n', total_water_used/num_zones);
fprintf('\nPer Zone Statistics:\n');

for zone = 1:num_zones
    zone_water = sum(water_consumed(zone, :));
    zone_events = sum(pump_status(zone, :));
    avg_moisture = mean(soil_moisture(zone, :));
    min_moisture = min(soil_moisture(zone, :));
    max_moisture = max(soil_moisture(zone, :));
    
    fprintf('Zone %d (%s):\n', zone, plant_profiles{zone, 1});
    fprintf('  Water used: %.2f L\n', zone_water);
    fprintf('  Irrigation events: %d\n', zone_events);
    fprintf('  Avg moisture: %.1f%%\n', avg_moisture);
    fprintf('  Min/Max moisture: %.1f%% / %.1f%%\n\n', min_moisture, max_moisture);
end

%% VISUALIZATION
% Figure 1: Soil Moisture Over Time
figure('Name', 'Soil Moisture Dynamics', 'Position', [100 100 1200 600]);
for zone = 1:num_zones
    subplot(2, 2, zone);
    
    % Plot moisture
    yyaxis left
    plot(time/24, soil_moisture(zone, :), 'b-', 'LineWidth', 2);
    hold on;
    
    % Plot threshold lines
    profile = plant_profiles{zone, 2};
    yline(profile(1), 'r--', 'Min Threshold', 'LineWidth', 1.5);
    yline(profile(2), 'g--', 'Max Threshold', 'LineWidth', 1.5);
    
    ylabel('Soil Moisture (%)');
    ylim([0 100]);
    
    % Plot irrigation events
    yyaxis right
    irrigation_times = time(pump_status(zone, :) == 1) / 24;
    stem(irrigation_times, ones(size(irrigation_times))*100, 'r', 'LineWidth', 1);
    ylabel('Irrigation Event');
    ylim([0 120]);
    
    xlabel('Time (days)');
    title(sprintf('Zone %d: %s', zone, plant_profiles{zone, 1}));
    grid on;
    legend('Soil Moisture', 'Min Threshold', 'Max Threshold', 'Irrigation');
end
sgtitle('Multi-Zone Soil Moisture and Irrigation Events', 'FontSize', 14, 'FontWeight', 'bold');

% Figure 2: Environmental Conditions
figure('Name', 'Environmental Conditions', 'Position', [150 150 1000 400]);
subplot(1, 2, 1);
plot(time/24, temperature, 'r-', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Temperature (°C)');
title('Ambient Temperature');
grid on;

subplot(1, 2, 2);
plot(time/24, humidity, 'b-', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Humidity (%)');
title('Ambient Humidity');
grid on;

% Figure 3: Water Consumption Analysis
figure('Name', 'Water Consumption', 'Position', [200 200 1000 500]);

% Cumulative water use
subplot(2, 2, 1);
cumulative_water = cumsum(sum(water_consumed, 1));
plot(time/24, cumulative_water, 'b-', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Cumulative Water (L)');
title('Total Water Consumption');
grid on;

% Water per zone (bar chart)
subplot(2, 2, 2);
zone_totals = sum(water_consumed, 2);
bar(zone_totals);
xlabel('Zone Number');
ylabel('Total Water (L)');
title('Water Consumption by Zone');
xticklabels(cellfun(@(x) x{1}, plant_profiles, 'UniformOutput', false));
xtickangle(45);
grid on;

% Irrigation frequency
subplot(2, 2, 3);
zone_events = sum(pump_status, 2);
bar(zone_events);
xlabel('Zone Number');
ylabel('Number of Events');
title('Irrigation Events by Zone');
xticklabels(cellfun(@(x) x{1}, plant_profiles, 'UniformOutput', false));
xtickangle(45);
grid on;

% Moisture control quality
subplot(2, 2, 4);
moisture_stats = [mean(soil_moisture, 2), std(soil_moisture, 0, 2)];
bar(moisture_stats);
xlabel('Zone Number');
ylabel('Moisture (%)');
title('Average Moisture ± Std Dev');
legend('Mean', 'Std Dev');
xticklabels(cellfun(@(x) x{1}, plant_profiles, 'UniformOutput', false));
xtickangle(45);
grid on;

%% SAVE RESULTS
results = struct();
results.time = time;
results.soil_moisture = soil_moisture;
results.temperature = temperature;
results.humidity = humidity;
results.pump_status = pump_status;
results.water_consumed = water_consumed;
results.plant_profiles = plant_profiles;

save('simulation_results.mat', 'results');
fprintf('\nResults saved to simulation_results.mat\n');

fprintf('\n=== SIMULATION COMPLETE ===\n');
