clc;
clear;


%==============================================
% mission timeline / deltaV maneuvers
%==============================================


deltaV1 = 15.00;
deltaV2 = 2162.55;
deltaV3 = 94.84;
deltaV4 = 237.3;
deltaV5 = 68.25;
deltaV6 = 20.00;
deltaV7 = 31.05;
deltaV8 = 10.00;
deltaV9 = 73.38;


sum_deltaV = deltaV1 + deltaV2 + deltaV3 + deltaV4 + deltaV5 + deltaV6 + deltaV7 + deltaV8 + deltaV9; % [m/s]

%=============================================
% Engines
%=============================================

%Thruster S10-26
p_tank_RCS = 20; % [bar]
F_RCS = ((12.5-6)/(23-10)) * p_tank_RCS + 1; % [N] thrust depending on tank pressure
m_engine_RCS = 16 * 0.65; % [kg] redundancy
ISP_RCS = 292; % [s]
OF_RCS = 1.65;
m_dot_RCS = (((4.2-2.3)/(23-10)) * p_tank_RCS + 0.838462) / 1000; % [kg/s] mass flow depending on pressure

%LEROS 4
F_bi = 1000; % [N]
m_engine_bi = 2 *8.41; % [kg]
ISP_bi = (321 * F_bi + 2* ISP_RCS * F_RCS) / (F_bi + 2 * F_RCS);  % [s] mixed ISP
OF_bi = 1.65;
p_tank_bi = p_tank_RCS*10^5; % [N/m^2]



%==============================================
% Loopstarter
%==============================================
m_wet = 1; % [kg]
m_wet_1 = 2; % [kg]
m_wet_2 = 1; % [kg]
m_prop = 1; % [kg]
t_1 = 1;
t_2 = 1;
t_3 = 1;
t_4 = 1;
t_5 = 1;
t_6 = 1;
t_7 = 1;
t_8 = 1;
t_9 = 1;

%==============================================
% power budget
%==============================================
p_PWR_Earth = 550; % [W]
p_PWR_Mars = p_PWR_Earth*1.524^2; % [W]


while abs(m_wet_1 - m_wet_2) > 10^-3

    %==============================================
    % propulsion budget
    %==============================================
    ISP_1 = ISP_bi; % [s] R-4D-15 HiPAT™ 445 N
    ISP_2 = ISP_bi; % [s]
    ISP_3 = ISP_bi; % [s]
    ISP_4 = ISP_bi; % [s]
    ISP_5 = ISP_bi; % [s]
    ISP_6 = ISP_bi; % [s]
    ISP_7 = ISP_bi; % [s]
    ISP_8 = ISP_bi; % [s]
    ISP_9 = ISP_bi; % [s]

   
    m_1_RCS = 2* m_dot_RCS * t_1; % [kg]
    m_2_RCS = 2* m_dot_RCS * t_2; % [kg]
    m_3_RCS = 2* m_dot_RCS * t_3; % [kg]
    m_4_RCS = 2* m_dot_RCS * t_4; % [kg]
    m_5_RCS = 2* m_dot_RCS * t_5; % [kg]
    m_6_RCS = 2* m_dot_RCS * t_6; % [kg]
    m_7_RCS = 2* m_dot_RCS * t_7; % [kg]
    m_8_RCS = 2* m_dot_RCS * t_8; % [kg]
    m_9_RCS = 2* m_dot_RCS * t_9; % [kg]
    m_prop_RCS = (m_1_RCS + m_2_RCS + m_3_RCS + m_4_RCS + m_5_RCS + m_6_RCS + m_7_RCS + m_8_RCS + m_9_RCS); % [kg] with margine
    m_MMH_RCS = 1/(OF_RCS+1) * m_prop_RCS;
    m_NTO_RCS = m_prop_RCS - m_MMH_RCS;

    roh_MMH = 880; % [kg/m^3]
    roh_NTO = 1440; % [kg/m^3]
    m_MMH = 1/(OF_bi+1) * m_prop; % [kg]
    m_NTO = m_prop - m_MMH; % [kg]
    V_MMH = (m_MMH / roh_MMH)*(1+0.10); % [m^3] pressure regulated
    V_NTO = (m_NTO / roh_NTO)*(1+0.10); % [m^3] pressure regulated
    r_i_MMH = (V_MMH / ((4/3)*pi)) ^ (1/3); % [m]
    r_i_NTO = (V_NTO / ((4/3)*pi)) ^ (1/3); % [m]
    sigma_material = 498*10^6; % [N/m^2] Titan
    roh_material = 4450; % [kg/m^3] Titan
    m_tank_MMH = 2*pi*r_i_MMH^3*p_tank_bi*roh_material/sigma_material * (1+0.993); % [kg]
    m_tank_NTO = 2*pi*r_i_NTO^3*p_tank_bi*roh_material/sigma_material * (1+0.993); % [kg]

    R_He = 2077.1; % [J/(kg*K)]
    T_He = 323.15; % [K]
    m_He = (p_tank_bi) * V_MMH / (R_He *T_He) + (p_tank_bi) * V_NTO / (R_He *T_He); % [kg]
    p_tank_He = 310*10^5; % [Pa] = [N/m^2]
    roh_He = p_tank_He / (R_He * T_He); % [kg/m^3] "He density at 310 bar at 20°C" Xenon Tank – COTS 
    V_He = (m_He / roh_He)*(1+0.20); % [m^3] 
    r_i_He = (V_He / ((4/3)*pi)) ^ (1/3); % [m]
    sigma_material_He = 733*10^6; % [N/m^2] CFRP
    roh_material_He = 1800; % [kg/m^3] CFRP
    m_tank_He = 2*pi*r_i_He^3*p_tank_He*roh_material_He/sigma_material_He * (1+0.752); % [kg] for Titanium liner


    m_PROP = (m_engine_bi + m_engine_RCS + m_tank_MMH + m_tank_NTO + m_tank_He) * (1+0.171); % [kg] with 17,1% margine for valves and so on
    
    %==============================================
    % mass budget
    %==============================================
    m_INS = 118.41; % [kg]
    m_AOGNC = 60.40; % [kg]
    m_COM = 71.64; % [kg]
    m_CDH = 38.48; % [kg]
    m_MEC = 39.00; % [kg]
    m_PWR = 350.04; % [kg]
    m_TC = 65.89; % [kg]
    m_STR = 303.26; % [kg]
    m_SHI = 237.02;
    m_RIN = 1.49;
    m_dry_no_Harness = m_INS + m_AOGNC + m_COM + m_CDH + m_MEC + m_PWR + m_TC + m_STR + m_SHI + m_RIN + m_PROP; % [kg]
    m_Harness = m_dry_no_Harness * 0.05; % [kg]
    m_dry = (m_dry_no_Harness + m_Harness) * (1+0.2); % [kg] with 20% margine 
    m_pres = m_He; % [kg] ullage volume
    m_wet_1 = m_dry + m_prop; % [kg]
    m_wet = m_wet_1 + m_pres; % [kg]
    m_adapter = 85; % [kg]
    m_launch = m_wet + m_adapter; % [kg]

    %==============================================
    % propellant budget
    %==============================================
    m_final_9 = m_dry; % [kg]
    m_0_9 = m_final_9 * exp((deltaV9/(9.81*ISP_9))); % [kg]
    m_prop_9 = m_0_9 - m_final_9; % [kg]
    m_dot_prop9 = F_bi/(ISP_9 * 9.81); % [kg/s]
    t_9 = (m_prop_9/m_dot_prop9); %[days]
    
    m_final_8 = m_0_9; % [kg]
    m_0_8 = m_final_8 * exp((deltaV8/(9.81*ISP_8))); % [kg]
    m_prop_8 = m_0_8 - m_final_8; % [kg]
    m_dot_prop8 = F_bi/(ISP_8 * 9.81); % [kg/s]
    t_8 = (m_prop_8/m_dot_prop8); %[days]

    m_final_7 = m_0_8; % [kg]
    m_0_7 = m_final_7 * exp((deltaV7/(9.81*ISP_7))); % [kg]
    m_prop_7 = m_0_7 - m_final_7; % [kg]
    m_dot_prop7 = F_bi/(ISP_7 * 9.81); % [kg/s]
    t_7 = (m_prop_7/m_dot_prop7); %[days]

    m_final_6 = m_0_7; % [kg]
    m_0_6 = m_final_6 * exp((deltaV6/(9.81*ISP_6))); % [kg]
    m_prop_6 = m_0_6 - m_final_6; % [kg]
    m_dot_prop6 = F_bi/(ISP_6 * 9.81); % [kg/s]
    t_6 = (m_prop_6/m_dot_prop6); %[days]

    m_final_5 = m_0_6; % [kg]
    m_0_5 = m_final_5 * exp((deltaV5/(9.81*ISP_5))); % [kg]
    m_prop_5 = m_0_5 - m_final_5; % [kg]
    m_dot_prop5 = F_bi/(ISP_5 * 9.81); % [kg/s]
    t_5 = (m_prop_5/m_dot_prop5); %[days]

    m_final_4 = m_0_5; % [kg]
    m_0_4 = m_final_4 * exp((deltaV4/(9.81*ISP_4))); % [kg]
    m_prop_4 = m_0_4 - m_final_4; % [kg]
    m_dot_prop4 = F_bi/(ISP_4 * 9.81); % [kg/s]
    t_4 = (m_prop_4/m_dot_prop4); %[days]

    m_final_3 = m_0_4; % [kg]
    m_0_3 = m_final_3 * exp((deltaV3/(9.81*ISP_3))); % [kg]
    m_prop_3 = m_0_3 - m_final_3; % [kg]
    m_dot_prop3 = F_bi/(ISP_3 * 9.81); % [kg/s]
    t_3 = (m_prop_3/m_dot_prop3); %[days]

    m_final_2 = m_0_3; % [kg]
    m_0_2 = m_final_2 * exp((deltaV2/(9.81*ISP_2))); % [kg]
    m_prop_2 = m_0_2 - m_final_2; % [kg]
    m_dot_prop2 = F_bi/(ISP_2 * 9.81); % [kg/s]
    t_2 = (m_prop_2/m_dot_prop2); %[days]

    m_final_1 = m_0_2; % [kg]
    m_0_1 = m_final_1 * exp((deltaV1/(9.81*ISP_1))); % [kg]
    m_prop_1 = m_0_1 - m_final_1; % [kg]
    m_dot_prop1 = F_bi/(ISP_1 * 9.81); % [kg/s]
    t_1 = (m_prop_1/m_dot_prop1); %[days]

    m_wet_2 = m_0_1; % [kg]
    m_prop = m_prop_9 +  m_prop_8 + m_prop_7 + m_prop_6 + m_prop_5 + m_prop_4 + m_prop_3 + m_prop_2 + m_prop_1; % [kg]
    

end

disp(["Wet mass [kg]",m_wet, "= 4398.33 kg"]);
disp(["Error [%]", abs(((4398.33/m_wet)-1)*100)]);
disp(["PROP Error [%]", abs(((233.94/m_PROP)-1)*100)]);
disp(["prop Error [%]", abs(((2469.06/m_prop)-1)*100)]);
disp(["Error without prop Error [%]", abs(((4398.33)/(m_wet-m_prop * abs(((2469.06/m_prop)-1)))-1)*100) ]);
disp(["D_V1 to D_V9 [m/s]",deltaV1, deltaV2, deltaV3, deltaV4, deltaV5, deltaV6, deltaV7, deltaV8, deltaV9, "Summ_Delta_V [m/s]", sum_deltaV]);
disp(["Dry mass [kg]", m_dry]);
disp(["PROP system mass [kg]", m_PROP]);
disp(["MMH tank mass [kg]", m_tank_MMH, "NTO tank mass [kg]", m_tank_NTO, "He tank mass [kg]", m_tank_He]);
disp(["MMH volume [l]", V_MMH*1000, "NTO volume [l]", V_NTO*1000, "He volume [l]", V_He*1000]);
disp(["MMH mass [kg]", m_MMH, "NTO mass [kg]", m_NTO, "He mass [kg]", m_He]);
disp(["Burn time[s]", t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9]);
disp(["m_0 [kg]",m_0_1, m_0_2, m_0_3, m_0_4, m_0_5, m_0_6, m_0_7, m_0_8, m_0_9]);
disp(["RCS prop mass [kg]", m_prop_RCS]);