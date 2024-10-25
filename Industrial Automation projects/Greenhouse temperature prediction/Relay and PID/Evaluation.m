%% 4th Point. Quadratic deviation for PID
PIDActualTemperature(:);
PIDDesiredTemperature = squeeze(PIDDesiredTemperature);
PIDDesiredTemperature = reshape(PIDDesiredTemperature, [], 1);
PIDdeviations = PIDActualTemperature - PIDDesiredTemperature;
PIDsquaredDeviations = PIDdeviations .^ 2;
PIDsumSquaredDeviations = sum(PIDsquaredDeviations);
PIDaverageSquaredDeviation = PIDsumSquaredDeviations / numel(PIDActualTemperature);
PIDquadraticDeviation = sqrt(PIDaverageSquaredDeviation);


% Quadratic deviation for Relay
RelayActualTemperature(:);
RelayDesiredTemperature = squeeze(RelayDesiredTemperature);
RelayDesiredTemperature = reshape(RelayDesiredTemperature, [], 1);
Relaydeviations = RelayActualTemperature - RelayDesiredTemperature;
RelaysquaredDeviations = Relaydeviations .^ 2;
RelaysumSquaredDeviations = sum(RelaysquaredDeviations);
RelayaverageSquaredDeviation = RelaysumSquaredDeviations / numel(RelayActualTemperature);
RelayquadraticDeviation = sqrt(RelayaverageSquaredDeviation);

disp("PIDquadraticDeviation");
disp(PIDquadraticDeviation);
disp("RelayquadraticDeviation");
disp(RelayquadraticDeviation);