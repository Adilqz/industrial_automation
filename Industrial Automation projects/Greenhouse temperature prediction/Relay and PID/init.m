% To remove all from the Workspace
clear; 
% To remove all from the Command Window
clc;
% Now we can start...
%==========================================
%% Module 1)
% Horizon 24 h = 24*60 min = 24*60*60 s
H = 24 * 60 * 60;% [s]

% Heating/cooling power range  
q_MAX_1 = 25000;% [W]
q_MAX = 20000;% [W]
q_min = 0;% [W]

% Below are the forecasted values for 24 hours for temperature and sun
% in Genoa for the 17th of April
% Define the solar radiation data in kWh/m2/d
solar_radiation = [0, 0, 0, 0, 0, 0, 2, 83, 335, 519, 657, ...
    751, 815, 809, 709, 632, 459, 208, 119, 1, 0, 0, 0, 0, 0];
solRad = solar_radiation(:);

% Define the temperature data in Kelvin
outside_temperature = [285, 289, 288, 288, 287, 286, 286, 285, 286, 288, 289, ...
    290, 291, 291, 292, 292, 291, 292, 291, 289, 289, 288, 287, 287, 286];
outTemp = outside_temperature(:);

desired_temperature = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 25, 25, 30, ...
    30, 30, 30, 30, 24, 23, 22, 22, 20, 20, 20, 20];
desired_temperature = Celsius2Kelvin(desired_temperature);
desTemp = desired_temperature(:);

% Define 24 hours
time = (0:1:24);
time = time(:);

% mass of the air in the module (let's consider that there is a 1.2 kg of 
% air per 1 cubic meter)
m_f = 1.2*300;% [Kg]
c_f = 1000;% [J/(Kg*K)] 
% Average thermal capacity of the mass of the air in the box
C_f = c_f*m_f;% [J/K]

% Average thermal transmittance per unit surface air/glass
barK_fc = 5.7; % [W/(m^2*K)]
% Surface of the glass
s_c = 100;% [m^2]
% Overall thermal transmittance between the air in the box and the glass
k_fc = barK_fc*s_c; % [W/K]

% Sun radiation
% the overall effective absorbance of the greenhouse
a = 0.97;
% The heat due to solar radiation transmitted through the greenhouse 
% covering material
q_s = barK_fc * a;
Q_s = barK_fc * a * solRad;

% Addictive "module" noise variance
var_wf = 1;

% 'Fromworkspace' block in Simulink needs timeseries 
%The '.*' operator is used to multiply each entry of the vector by 60*60
%(indeed it's given in hours, but SI time unit is second s). Transpose
%since we want a row vector, while it's a column vector.
timeseriesOutTemp = timeseries(outTemp', time.*(60*60)');
timeseriesSolRad = timeseries(solRad', time.*(60*60)');
timeseriesDesTemp = timeseries(desTemp', time.*(60*60)');

% Sample time
Ts = 5;%[s]

theta0_f = 280;

%% Module 2)

% mass of the air in the module (let's consider that there is a 1.2 kg of 
% air per 1 cubic meter)
m_f_2 = 1.2*300;% [Kg]
c_f_2 = 1000;% [J/(Kg*K)] 
% Average thermal capacity of the mass of the air in the box
C_f_2 = c_f_2*m_f_2;% [J/K]

%% Module 3)

% mass of the air in the module (let's consider that there is a 1.2 kg of 
% air per 1 cubic meter)
m_f_3 = 1.2*200;% [Kg]
c_f_3 = 1000;% [J/(Kg*K)] 
% Average thermal capacity of the mass of the air in the box
C_f_3 = c_f_3*m_f_3;% [J/K]

%% Module 4)

% mass of the air in the module (let's consider that there is a 1.2 kg of 
% air per 1 cubic meter)
m_f_4 = 1.2*400;% [Kg]
c_f_4 = 1000;% [J/(Kg*K)] 
% Average thermal capacity of the mass of the air in the box
C_f_4 = c_f_4*m_f_4;% [J/K]

%% 
q_star = 5000;




%% Power consumption


%% Funtions

function kelvin = Celsius2Kelvin(celsius)
    kelvin = celsius + 273;
end






