clc;
clear;


m_dry = 1121.5894 + 3.5251 + 20; % [kg] dry + helium + p_probe
m_probe = 20; % [kg]

% ------ ELECTRIC ---------
t_9 = 9035872.43670; % [s] EOL
dV_8 = 104.9497; % [m/s] Maintenance
% arrive at Phobos Orbit
t_7 = 11334710.1876; % [s] Spiral Down
t_6 = 10164236.2100; % [s] Circularize


% ------ CHEMICAL ---------
t_5 =   13.1293586301;  % [s] Raise Peri
t_4 =   42.4428661248;  % [s] MatchPlane
t_3 = 2858.30549682;    % [s] MOI
t_2 =    5.34999976365; % [s] TCM
t_1 = 2155.84586111;    % [s] C3


ISP_bi = 319.7364; % [s]
ISP_el = 1800; % [s]
F_bi = 449.9997; % [N]
F_el = 0.090; % [N]
v_exit_bi = ISP_bi * 9.81; % [m/s]
v_exit_el = ISP_el * 9.81; % [m/s]

m_dot_bi = F_bi/(ISP_bi * 9.81); % [kg/s]
m_dot_el = F_el/(ISP_el * 9.81); % [kg/s]

m_prop_1 = m_dot_bi * t_1; % [kg]
m_prop_2 = m_dot_bi * t_2; % [kg]
m_prop_3 = m_dot_bi * t_3; % [kg]
m_prop_4 = m_dot_bi * t_4; % [kg]
m_prop_5 = m_dot_bi * t_5; % [kg]
m_prop_6 = m_dot_el * t_6; % [kg]
m_prop_7 = m_dot_el * t_7; % [kg]
m_prop_9 = m_dot_el * t_9; % [kg]

m_0_9 = m_dry + m_prop_9 - m_probe; % [kg]

m_0_8 = m_0_9 * exp((dV_8/(9.81*ISP_el))); % [kg]
m_prop_8 = m_0_8 - m_0_9; % [kg]

m_0_7 = m_dry + m_prop_9 + m_prop_8 + m_prop_7; % [kg]
m_0_6 = m_dry + m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6; % [kg]
m_0_5 = m_dry + m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6 + m_prop_5; % [kg]
m_0_4 = m_dry + m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6 + m_prop_5 + m_prop_4; % [kg]
m_0_3 = m_dry + m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6 + m_prop_5 + m_prop_4 + m_prop_3; % [kg]
m_0_2 = m_dry + m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6 + m_prop_5 + m_prop_4 + m_prop_3 + m_prop_2; % [kg]
m_0_1 = m_dry + m_prop_9 + m_prop_8 + m_prop_7 + m_prop_6 + m_prop_5 + m_prop_4 + m_prop_3 + m_prop_2 + m_prop_1; % [kg]

dV_9 =  log(((-t_9*F_el)/(m_0_9*v_exit_el))+1) * v_exit_el * (-1); % [m/s]
dV_7 =  log(((-t_7*F_el)/(m_0_7*v_exit_el))+1) * v_exit_el * (-1); % [m/s]
dV_6 =  log(((-t_6*F_el)/(m_0_6*v_exit_el))+1) * v_exit_el * (-1); % [m/s]
dV_5 =  log(((-t_5*F_bi)/(m_0_5*v_exit_bi))+1) * v_exit_bi * (-1); % [m/s]
dV_4 =  log(((-t_4*F_bi)/(m_0_4*v_exit_bi))+1) * v_exit_bi * (-1); % [m/s]
dV_3 =  log(((-t_3*F_bi)/(m_0_3*v_exit_bi))+1) * v_exit_bi * (-1); % [m/s]
dV_2 =  log(((-t_2*F_bi)/(m_0_2*v_exit_bi))+1) * v_exit_bi * (-1); % [m/s]
dV_1 =  log(((-t_1*F_bi)/(m_0_1*v_exit_bi))+1) * v_exit_bi * (-1); % [m/s]

dV_chem = dV_1 + dV_2 + dV_3 + dV_4 + dV_5;
dV_el = dV_6 + dV_7 + dV_8 + dV_9;
dV_tot = dV_chem + dV_el;

m_prop_chem = m_prop_1 + m_prop_2 + m_prop_3 + m_prop_4 + m_prop_5;
m_prop_el = m_prop_6 + m_prop_7 + m_prop_8 + m_prop_9;


disp(["dV [m/s]", dV_1, dV_2, dV_3, dV_4, dV_5, dV_6 ,dV_7, dV_8, dV_9]);
disp(["dV tot [m/s]", dV_tot, "dV chem [m/s]", dV_chem, "dV el [m/s]", dV_el]);
disp(["m_prop_1 [kg]",m_prop_1]);
disp(["m_prop_2 [kg]",m_prop_2]);
disp(["m_prop_3 [kg]",m_prop_3]);
disp(["m_prop_4 [kg]",m_prop_4]);
disp(["m_prop_5 [kg]",m_prop_5]);
disp(["m_prop_6 [kg]",m_prop_6]);
disp(["m_prop_7 [kg]",m_prop_7]);
disp(["m_prop_8 [kg]",m_prop_8]);
disp(["m_prop_9 [kg]",m_prop_9]);
disp(["m_prop_chem [kg]",m_prop_chem]);
disp(["m_prop_el [kg]",m_prop_el]);

