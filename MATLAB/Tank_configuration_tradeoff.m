clc;

%==============================================
% Given volumes of pressurant and propellants
%==============================================

%Configuration tradeoff considering approximate total needs of propellant
%and pressurant:
% V_He= 46.6 l
% V_MMH_GMAT= 322.2 l
% V_MON_GMAT= 325.3l
% V_Noble_GMAT = 58.6 l

%==============================================
% Available tank mass and volume characteristics
%==============================================

% Propellant tank 1 Ariane 198 litre, 282 litre, 331 litre bipropellant surface tension tank. Model OST 25-0
V_Proptank_Ariane25(1) = 0.198; %[m3]
V_Proptank_Ariane25(2) = 0.282; %[m3]
V_Proptank_Ariane25(3) = 0.331; %[m3]
m_Proptank_Ariane25(1) = 21+((36-21)/(700-282))*((V_Proptank_Ariane25(1)*1000)-282);
m_Proptank_Ariane25(2) = 21+((36-21)/(700-282))*((V_Proptank_Ariane25(2)*1000)-282); %[kg]
m_Proptank_Ariane25(3) = 21+((36-21)/(700-282))*((V_Proptank_Ariane25(3)*1000)-282); %[kg]
p_Proptank_Ariane25 = 33e5; %[Pa]
T_Proptank_Ariane25= 293.15; %[K]

% Propellant tank 2 235 Litre Bipropellant Tank. Model OST 21/0
V_Proptank_Ariane21= 0.235; %[m3]
m_Proptank_Ariane21= 16; %[kg]
p_Proptank_Ariane21= 22e5; % [Pa]
T_Proptank_Ariane21= 293.15; %[K]


%Propellant tank 3 700 to 1108 Litre Bipropellant Tank. Model OST 22/X
V_Proptank_Ariane22(1)=0.7; %[m3]
m_Proptank_Ariane22(1) = 36+((49-36)/(1108-700))*((V_Proptank_Ariane22(1)*1000)-700); %[kg]
p_Proptank_Ariane22 = 19.5e5;
T_Proptank_Ariane22= 293.15; %[K]


% Helium tank 1 PVG Family 40–75L High Pressure Tank (HPV) – Helium/Nitrogen – COTS
V_HeTank(1)= 0.040; % [m3]
m_HeTank(1) = 8.5+((14.4-8.5)/(75-40))*((V_HeTank(1)*1000)-40); %[kg]
V_HeTank(2)= 0.070; %[m3]
m_HeTank(2) = 8.5+((14.4-8.5)/(75-40))*((V_HeTank(2)*1000)-40); %[kg]

%Helium tank 2 PVG Family 80-120 High Pressure Tank (HPV) – Helium/Nitrogen – COTS with delta Qual
V_HeTank(3)= 0.095; %[kg] Any value between 80 and 120
m_HeTank(3) = 8.5+((14.4-8.5)/(75-40))*((V_HeTank(2)*1000)-40); %[kg]
%m_HeTank(2)=23.5; %[kg]
p_HeTank =310e5; %[Pa]
T_HeTank = 393.15; %[K]

%==============================================
% Configuration 1: 2 tanks MMH, 2 tanks MON, 1 tank He
%==============================================

%Initial free volumes in propellant tanks
V_MMHTank_Ullage(1) = 2*V_Proptank_Ariane25(1) - V_MMH_GMAT; 
V_MONTank_Ullage(1) = 2*V_Proptank_Ariane25(1) - V_MON_GMAT; 
%Density of helium in propellant tanks
roh_He_MMH(1) = py.CoolProp.CoolProp.PropsSI('D','P',p_Proptank_Ariane25,'T',T_Proptank_Ariane25,'Helium'); % [kg/m^3]
roh_He_MON(1) = py.CoolProp.CoolProp.PropsSI('D','P',p_Proptank_Ariane25,'T',T_Proptank_Ariane25,'Helium'); % [kg/m^3]
%Initial additional helium mass to fill up propellant tank ullage volume
m_He_initial(1) = V_MMHTank_Ullage(1)*roh_He_MMH(1)+V_MONTank_Ullage(1)*roh_He_MON(1); % [kg]

%Total dry mass of configuration
m_config_dry_MON(1)= m_HeTank(2)+2* m_Proptank_Ariane25(1) +2*m_Proptank_Ariane25(1);
%Total mass of configuration considering only additional Helium and Xenon
m_config_wet_MON(1)= m_config_dry_MON(1)+ m_He_initial(1);

%==============================================
% Configuration 2: 2 tanks MMH, 2 tanks MON, 2 tank He
%==============================================

%Initial free volumes in propellant tanks
V_MMHTank_Ullage(2) = V_MMHTank_Ullage(1); 
V_MONTank_Ullage(2) = V_MONTank_Ullage(1); 
%Density of helium in propellant tanks
roh_He_MMH(2) = roh_He_MMH(1); % [kg/m^3]
roh_He_MON(2) = roh_He_MON(1); % [kg/m^3]
%Initial additional helium mass to fill up propellant tank ullage volume
m_He_initial(2) = m_He_initial(1); % [kg]
%Volume of Helium in pressurant tank
roh_He_He(2) = roh_He_He(1); % [kg/m^3]


%Total dry mass of configuration
m_config_dry_MON(2)= 2*m_HeTank(1)+2* m_Proptank_Ariane25(1) +2*m_Proptank_Ariane25(1);
%Total mass of configuration considering only additional Helium and extra
%Helium
m_config_wet_MON(2)= m_config_dry_MON(2)+ m_He_initial(1);

%==============================================
% Configuration 3: 1 tanks MMH, 2 tanks MON, 2 tanks He
%==============================================
%Initial free volumes in propellant tanks
V_MMHTank_Ullage(3) = V_Proptank_Ariane22(1) - V_MMH_GMAT; 
V_MONTank_Ullage(3) = 2*V_Proptank_Ariane25(1) - V_MON_GMAT; 
%Density of helium in propellant tanks
roh_He_MMH(3) = py.CoolProp.CoolProp.PropsSI('D','P',p_Proptank_Ariane25,'T',T_Proptank_Ariane25,'Helium'); % [kg/m^3]
roh_He_MON(3) = py.CoolProp.CoolProp.PropsSI('D','P',p_Proptank_Ariane25,'T',T_Proptank_Ariane25,'Helium'); % [kg/m^3]
%Initial additional helium mass to fill up propellant tank ullage volume
m_He_initial(3) = V_MMHTank_Ullage(3)*roh_He_MMH(3)+V_MONTank_Ullage(3)*roh_He_MON(3); % [kg]
%Volume of Helium in pressurant tank
roh_He_He(3) = py.CoolProp.CoolProp.PropsSI('D','P',p_HeTank,'T',T_HeTank,'Helium'); % [kg/m^3]

%Total dry mass of configuration
m_config_dry_MON(3)= 2*m_HeTank(1)+m_Proptank_Ariane25(3) +2*m_Proptank_Ariane25(1);
%Total mass of configuration considering only additional Helium
m_config_wet_MON(3)= m_config_dry_MON(3)+ m_He_initial(3);

%==============================================
% Configuration 4: 1 tanks MMH, 1 tanks MON, 1 tanks He
%==============================================
%Initial free volumes in propellant tanks
V_MMHTank_Ullage(4) = V_Proptank_Ariane25(3) - V_MMH_GMAT; 
V_MONTank_Ullage(4) = V_Proptank_Ariane22(1) - V_MON_GMAT; 
%Density of helium in propellant tanks
roh_He_MMH(4) = py.CoolProp.CoolProp.PropsSI('D','P',p_Proptank_Ariane25,'T',T_Proptank_Ariane25,'Helium'); % [kg/m^3]
roh_He_MON(4) = py.CoolProp.CoolProp.PropsSI('D','P',p_Proptank_Ariane22,'T',T_Proptank_Ariane22,'Helium'); % [kg/m^3]
%Initial additional helium mass to fill up propellant tank ullage volume
m_He_initial(4) = V_MMHTank_Ullage(4)*roh_He_MMH(4)+V_MONTank_Ullage(4)*roh_He_MON(4); % [kg]
%Volume of Helium in pressurant tank
roh_He_He(4) = py.CoolProp.CoolProp.PropsSI('D','P',p_HeTank,'T',T_HeTank,'Helium'); % [kg/m^3]

%Total dry mass of configuration
m_config_dry_MON(4)= m_HeTank(2)+m_Proptank_Ariane25(3) +m_Proptank_Ariane22(1);
%Total mass of configuration considering only additional Helium
m_config_wet_MON(4)= m_config_dry_MON(4)+ m_He_initial(4);

disp(["Dry mass tanks ",m_config_dry_MON]);
disp(["Wet mass tanks ",m_config_wet_MON]);

% Create figure
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(m_config_dry_MON,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0.850980392156863 0.325490196078431 0.0980392156862745]);
plot(m_config_wet_MON,'MarkerSize',8,'Marker','+','LineWidth',2,'LineStyle','none','Color',[0 0 1]);
ylabel('Mass[kg]');
xlabel('Configuration option');
title('Configuration choice');
box(axes1,'on');
hold(axes1,'off');
set(axes1,'XTick',[1 2 3 4]);
