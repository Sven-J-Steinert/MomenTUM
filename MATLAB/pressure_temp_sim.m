clc;
%==============================================
% Pressurization Tank Sizing
%==============================================
%==============================================
% Read GMAT Values
%==============================================
% Load the data from the Excel file into a table
filename = 'GMAT_values_16.xlsx';
data_table= readmatrix(filename);

%Maneuver Time GMAT Values
t_man(1)= 0;
t_man(2)= data_table(1,3);
t_man(3)= data_table(2,4);
t_man(4)= data_table(3,5);
t_man(5)= data_table(4,6);
t_man(6)= data_table(5,7);

%MMH mass required every maneuver
m_MMH_maneuver(1)=0;
m_MMH_maneuver(2)=data_table(1,16)-data_table(2,16);
m_MMH_maneuver(3)=data_table(2,16)-data_table(3,16);
m_MMH_maneuver(4)=data_table(3,16)-data_table(4,16);
m_MMH_maneuver(5)=data_table(4,16)-data_table(5,16);
m_MMH_maneuver(6)=data_table(5,16)-data_table(6,16);

%MON mass required every maneuver
m_MON_maneuver(1)=0;
m_MON_maneuver(2)=data_table(1,17)-data_table(2,17);
m_MON_maneuver(3)=data_table(2,17)-data_table(3,17);
m_MON_maneuver(4)=data_table(3,17)-data_table(4,17);
m_MON_maneuver(5)=data_table(4,17)-data_table(6,17);
m_MON_maneuver(6)=data_table(5,17)-data_table(6,17);

%Propellant mass required every maneuver (Values from GMAT)
m_propellant(1) = 0;
m_propellant(2)= m_MMH_maneuver(2) + m_MON_maneuver(2);
m_propellant(3)= m_MMH_maneuver(3) + m_MON_maneuver(3);
m_propellant(4)= m_MMH_maneuver(4) + m_MON_maneuver(4);
m_propellant(5)= m_MMH_maneuver(5) + m_MON_maneuver(5);
m_propellant(6)= m_MMH_maneuver(6) + m_MON_maneuver(6);

m_propellant_total= m_propellant(1)+m_propellant(2)+m_propellant(3)+m_propellant(4)+m_propellant(5)+m_propellant(6);
disp(m_propellant_total);
disp(m_propellant)

%Initial MMH Mass and Volume GMAT
m_MMH_GMAT = data_table(1,16);
V_MMH_GMAT = (m_MMH_GMAT / roh_MMH); % [m^3] pressure regulated
%Initial MON Mass and Volume GMAT
m_MON_GMAT = data_table(1,17);
V_MON_GMAT = (m_MON_GMAT / roh_MON); % [m^3] pressure regulated

%==============================================
% MMH Tank
%==============================================
V_tank_MMH = 0.331; % [m^3] Tank 198l MEOP=20 bar T=293,15K
p_tank_MMH= 20e5; %[Pa]
T_tank_MMH = 293.15; %[K]
roh_He_MMHTank(1) = py.CoolProp.CoolProp.PropsSI('D','P',p_tank_MMH,'T',T_tank_MMH,'Helium'); % [kg/m^3]

%==============================================
% MON Tank
%==============================================
V_tank_MON = 0.198*2; % [m^3] Tank 198l MEOP=20 bar T=293,15K
p_tank_MON=20e5; %[Pa]
T_tank_MON = 293.15; %[K]
roh_He_MONTank(1) =  py.CoolProp.CoolProp.PropsSI('D','P',p_tank_MON,'T',T_tank_MON,'Helium'); % [kg/m^3]

%==============================================
% Propellant Values
%==============================================
%MMH Values from Sutton2010_Rocket_Propulsion_Elements_8th_Edition_by_Oscar_Biblarz_George_P._Sutton.pdf
p_vap_MMH= 0.638e5; %[Pa] at 67 degrees 
T_sat_MMH= 340.15; %[K] at 20 bar pressure

%MON-1 chosen due to Thruster requirement and to facilitate storage
%conditions (CEA Values)
p_vap_MON=0.0009e5; %[Pa] at -17 degrees
T_sat_MON= 296.15; %[K] at 20 bar pressure

%==============================================
% Pressurant Values
%==============================================
cp_He= 20.79; % (J/mol*K)
M_He= 4.002602e-3; %[kg/mol]
cp_He_kg=cp_He/M_He; %[J/kg*K]
JT_coeff =-0.060; %[K/bar]

%==============================================
% Vector Definition
%==============================================

%Total Maneuver Time 6 points
t_man_tot0(1)=0;
t_man_tot0(2)=t_man(2);
t_man_tot0(3)=t_man(2)+t_man(3);
t_man_tot0(4)=t_man(2)+t_man(3)+t_man(4);
t_man_tot0(5)=t_man(2)+t_man(3)+t_man(4)+t_man(5);
t_man_tot0(6)=t_man(2)+t_man(3)+t_man(4)+t_man(5)+t_man(6);

%Total Maneuver Time 9 points
t_wait=0;
t_man_tot(1)=0;
t_man_tot(2)=t_man(2);
t_man_tot(3)=t_man_tot(2)+t_man(3);
t_man_tot(4)=t_man_tot(3)+t_wait;
t_man_tot(5)=t_man_tot(4)+t_man(4);
t_man_tot(6)=t_man_tot(5)+t_wait;
t_man_tot(7)=t_man_tot(6)+t_man(5);
t_man_tot(8)=t_man_tot(7)+t_wait;
t_man_tot(9)=t_man_tot(8)+t_man(6);

%==============================================
% Required Helium in MMH Tank
%==============================================
%MMH Mass every maneuver
%m_MMH_maneuver = m_propellant/(OF_bi+1);

%MMH Volume every maneuver
%Assumption: MMH has constant p and T conditions
V_MMH_maneuver = m_MMH_maneuver/roh_MMH;
%MMH Volume flow every maneuver
Vdot_MMH = V_MMH_maneuver ./ t_man;
mdot_MMH= m_MMH_maneuver./ t_man;

%==============================================
% Required Helium in MON Tank
%==============================================
%MON Mass every maneuver 
%m_MON_maneuver = m_propellant*(OF_bi/(OF_bi+1));

%MON Volume every maneuver
%Assumption: MMH has constant p and T conditions
V_MON_maneuver = m_MON_maneuver/roh_MON;
%MON Volume flow every maneuver
Vdot_MON = V_MON_maneuver ./ t_man;
mdot_MON= m_MON_maneuver./ t_man;

%Helium Volume required each maneuver
V_He_Maneuver = V_MMH_maneuver + V_MON_maneuver;


%Mass of Helium required for first maneuver to fill volume in MMH Tank
%assuming both MMH and NTO tank have the same Temperature and Pressure
%conditions
%Mass of Helium for MMH Tank and MON tank for point 0
m_He_maneuver_MMHTank(1)= V_MMH_maneuver(1)* roh_He_MMHTank(1);
m_He_maneuver_MONTank(1)=V_MON_maneuver(1)* roh_He_MONTank(1);
m_He_maneuver_MMHTank_final1(1)=0;
m_He_maneuver_MMHTank_final2(1)=0;

%==============================================
% Helium Initial State in Helium Tank
%==============================================
P(1) = 310*10^5;
T(1)= 323.15;
T_const=323.15;
V_tank_He = 2*0.04; % [m^3]
roh_He_HeTank(1) = py.CoolProp.CoolProp.PropsSI('D','P',P(1),'T',T(1),'Helium'); % [kg/m^3]
m_He_HeTank(1) = roh_He_HeTank(1) * V_tank_He; % [kg]
H(1)= py.CoolProp.CoolProp.PropsSI('H','P',P(1),'T',T(1),'Helium'); % [J/kg] 
U(1) = py.CoolProp.CoolProp.PropsSI('U','P',P(1),'T',T(1),'Helium'); % [J/kg] 
U_He(1) = U(1)*m_He_HeTank(1);

%==============================================
% Helium Initial State in MMH Tank
%==============================================
p_He_MMHTank = p_tank_MMH;
%Ullage Volume
V_ullage_MMH = V_tank_MMH-V_MMH_GMAT;
V_He_MMHTank(1) = V_ullage_MMH;
T_He_MMHTank(1) = 293.15; %[K]
roh_He_MMHTank(1) = py.CoolProp.CoolProp.PropsSI('D','P',p_He_MMHTank,'T',T_He_MMHTank(1),'Helium'); % [kg/m^3]
%Initial mass of Helium in MMH Tank
m_He_MMHTank(1) = V_ullage_MMH*roh_He_MMHTank(1);

%==============================================
% Helium Initial State in MON Tank
%==============================================
p_He_MONTank = p_tank_MON;
%Ullage Volume
V_ullage_MON = V_tank_MON-V_MON_GMAT;
V_He_MONTank(1) = V_ullage_MON;
T_He_MONTank(1) = 293.15; %[K]
roh_He_MONTank(1) = py.CoolProp.CoolProp.PropsSI('D','P',p_He_MONTank,'T',T_He_MONTank(1),'Helium'); % [kg/m^3]
%Initial mass of Helium in MON Tank
m_He_MONTank(1) = V_ullage_MON*roh_He_MONTank(1);

%Initial conditions Helium Mass for Maneuver 1
m_He_maneuver_MMHTank(2)= V_MMH_maneuver(2)* roh_He_MMHTank(1);
m_He_maneuver_MONTank(2)=V_MON_maneuver(2)* roh_He_MONTank(1);

%==============================================
% Helium Pressure Drop Calculation
%==============================================
%Heat Flow during Maneuvers
Q_dot(1) = 50; % [W]
Q_dot(2) = 0; % [W]
Q_dot(3) = 250; % [W]
Q_dot(4) = 0; % [W]
Q_dot(5) = 0; % [W]


%==============================================
% Maneuver 1
%==============================================
disp("Maneuver 1");

%Helium in Propellant Tanks
%Temperature increase of Helium over Pressure regulator initial value
dT(1) = JT_coeff*(p_He_MMHTank*10^-5-P(1)*10^-5);
T_He_PR(1)=T(1)+dT(1);

for i=1:2

%MMH Tank
%Initial values
m_He_maneuver_MMHTank_final1(2)=m_He_maneuver_MMHTank(2);
m_He_maneuver_MMHTank_final2(2)=1;
while abs(m_He_maneuver_MMHTank_final1(2)-m_He_maneuver_MMHTank_final2(2))>10^-3
T_He_MMHTank(2)= (cp_He_kg*(m_He_MMHTank(1)*T_He_MMHTank(1)+m_He_maneuver_MMHTank_final1(2)*T_He_PR(1)))/(cp_He_kg*(m_He_MMHTank(1)+m_He_maneuver_MMHTank_final1(2)));
roh_He_MMHTank(2)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MMHTank,'T',T_He_MMHTank(2), 'Helium'); % [J/kg]
m_He_maneuver_MMHTank_final2(2)=V_MMH_maneuver(2)*roh_He_MMHTank(2);
m_He_maneuver_MMHTank_final1(2)=m_He_maneuver_MMHTank_final1(2)-0.001;
%disp(["Error", abs(m_He_maneuver_MMHTank_final1(2)-m_He_maneuver_MMHTank_final2(2))]);
end



%MON Tank
%Initial values
m_He_maneuver_MONTank_final1(2)=m_He_maneuver_MONTank(2);
m_He_maneuver_MONTank_final2(2)=1;
while abs(m_He_maneuver_MONTank_final1(2)-m_He_maneuver_MONTank_final2(2))>10^-3
T_He_MONTank(2)= (cp_He_kg*(m_He_MONTank(1)*T_He_MONTank(1)+m_He_maneuver_MONTank_final1(2)*T_He_PR(1)))/(cp_He_kg*(m_He_MONTank(1)+m_He_maneuver_MONTank_final1(2)));
roh_He_MONTank(2)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MONTank,'T',T_He_MONTank(2), 'Helium'); % [J/kg]
m_He_maneuver_MONTank_final2(2)=V_MON_maneuver(2)*roh_He_MONTank(2);
m_He_maneuver_MONTank_final1(2)=m_He_maneuver_MONTank_final1(2)-0.001;
%disp(["Error", abs(m_He_maneuver_MONTank_final1(2)-m_He_maneuver_MONTank_final2(2))]);
end

%Total Helium required for maneuver 1
m_He_maneuver_final(2)= m_He_maneuver_MONTank_final2(2)+m_He_maneuver_MMHTank_final2(2);
mdot_He_maneuver_final(2)= m_He_maneuver_final(2)/t_man(2);
m_He_HeTank (2)= m_He_HeTank (1) - m_He_maneuver_final(2);
roh_He_HeTank(2) = m_He_HeTank(2)/V_tank_He; 

%Helium in Helium Tank
U_He(2) = U_He(1)+t_man(2)*(Q_dot(1)-H(1)*mdot_He_maneuver_final(2)); % [kJ]
U(2) = U_He(2)/m_He_HeTank(2); % [kJ/kg]

P(2) = py.CoolProp.CoolProp.PropsSI('P','U',U(2),'D',roh_He_HeTank(2), 'Helium'); % [J/kg]
T(2) = py.CoolProp.CoolProp.PropsSI('T','U',U(2),'D',roh_He_HeTank(2), 'Helium'); % [J/kg]
H(2) = py.CoolProp.CoolProp.PropsSI('H','U',U(2),'D',roh_He_HeTank(2), 'Helium'); % [J/kg]



dT(1) = JT_coeff*(p_He_MMHTank*10^-5-((P(1)+P(2))/2)*10^-5);
T_He_PR(1)=(T(1)+T(2))/2+dT(1);
end

%Helium left in Propellant Tanks after first maneuver
m_He_MMHTank(2)= m_He_MMHTank(1)+m_He_maneuver_MMHTank_final2(2);
m_He_MONTank(2)=m_He_MONTank(1)+m_He_maneuver_MONTank_final2(2);

%Initial values in propellant tanks for second maneuver
m_He_maneuver_MMHTank(3)= V_MMH_maneuver(3)* roh_He_MMHTank(2);
m_He_maneuver_MONTank(3)=V_MON_maneuver(3)* roh_He_MONTank(2);

%==============================================
% Maneuver 2
%==============================================
disp("Maneuver 2");
%Helium in Propellant Tanks
%Temperature increase of Helium over Pressure regulator initial value
dT(2) = JT_coeff*(p_He_MMHTank*10^-5-P(2)*10^-5);
T_He_PR(2)=T(2)+dT(2);

for i=1:2
%MMH Tank
%Initial values
m_He_maneuver_MMHTank_final1(3)=m_He_maneuver_MMHTank(3);
m_He_maneuver_MMHTank_final2(3)=1;

while abs(m_He_maneuver_MMHTank_final1(3)-m_He_maneuver_MMHTank_final2(3))>10^-3
T_He_MMHTank(3)= (cp_He_kg*(m_He_MMHTank(2)*T_He_MMHTank(2)+m_He_maneuver_MMHTank_final1(3)*T_He_PR(2)))/(cp_He_kg*(m_He_MMHTank(2)+m_He_maneuver_MMHTank_final1(3)));
roh_He_MMHTank(3)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MMHTank,'T',T_He_MMHTank(3), 'Helium'); % [J/kg]
m_He_maneuver_MMHTank_final2(3)=V_MMH_maneuver(3)*roh_He_MMHTank(3);
m_He_maneuver_MMHTank_final1(3)=m_He_maneuver_MMHTank_final1(3)+0.001;
%disp(["Error", abs(m_He_maneuver_MMHTank_final1(3)-m_He_maneuver_MMHTank_final2(3))]);
end

%MON Tank
%Initial values
m_He_maneuver_MONTank_final1(3)=m_He_maneuver_MONTank(3);
m_He_maneuver_MONTank_final2(3)=1;
while abs(m_He_maneuver_MONTank_final1(3)-m_He_maneuver_MONTank_final2(3))>10^-3
T_He_MONTank(3)= (cp_He_kg*(m_He_MMHTank(2)*T_He_MONTank(2)+m_He_maneuver_MONTank_final1(3)*T_He_PR(2)))/(cp_He_kg*(m_He_MONTank(2)+m_He_maneuver_MONTank_final1(3)));
roh_He_MONTank(3)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MONTank,'T',T_He_MONTank(3), 'Helium'); % [J/kg]
m_He_maneuver_MONTank_final2(3)=V_MON_maneuver(3)*roh_He_MONTank(3);
m_He_maneuver_MONTank_final1(3)=m_He_maneuver_MONTank_final1(3)+0.001;
%disp(["Error", abs(m_He_maneuver_MONTank_final1(3)-m_He_maneuver_MONTank_final2(3))]);
end

%Total Helium required for maneuver 2
m_He_maneuver_final(3)= m_He_maneuver_MONTank_final2(3)+m_He_maneuver_MMHTank_final2(3);
mdot_He_maneuver_final(3)= m_He_maneuver_final(3)/t_man(3);
m_He_HeTank (3)= m_He_HeTank (2) - m_He_maneuver_final(3);
roh_He_HeTank(3) = m_He_HeTank(3)/V_tank_He ; 

%Helium in Helium Tank
U_He(3) = U_He(2)+t_man(3)*(Q_dot(2)-H(2)*mdot_He_maneuver_final(3)); % [kJ]
U(3) = U_He(3)/m_He_HeTank(3); % [kJ/kg]

P(3) = py.CoolProp.CoolProp.PropsSI('P','U',U(3),'D',roh_He_HeTank(3), 'Helium'); % [J/kg]
T(3) = py.CoolProp.CoolProp.PropsSI('T','U',U(3),'D',roh_He_HeTank(3), 'Helium'); % [J/kg]
H(3) = py.CoolProp.CoolProp.PropsSI('H','U',U(3),'D',roh_He_HeTank(3), 'Helium'); % [J/kg]


%Temperature increase of Helium over Pressure regulator
dT(2) = JT_coeff*(p_He_MMHTank*10^-5-((P(2)+P(3))/2)*10^-5);
T_He_PR(2)=(T(2)+T(3))/2+dT(2);
end

%Helium left in Propellant Tanks after second maneuver
m_He_MMHTank(3)= m_He_MMHTank(2)+m_He_maneuver_MMHTank_final2(3);
m_He_MONTank(3)=m_He_MONTank(2)+m_He_maneuver_MONTank_final2(3);

%Initial values in propellant tanks for third maneuver
m_He_maneuver_MMHTank(4)= V_MMH_maneuver(4)* roh_He_MMHTank(3);
m_He_maneuver_MONTank(4)=V_MON_maneuver(4)* roh_He_MONTank(3);

%==============================================
% Waiting time between Maneuver 2 and 3
%==============================================
T(4) = T_const; % [5]
P(4) = py.CoolProp.CoolProp.PropsSI('P','T',T(4),'D',roh_He_HeTank(3), 'Helium'); % [J/kg]
H(4)= py.CoolProp.CoolProp.PropsSI('H','P',P(4),'T',T(4), 'Helium'); % [J/kg]
U(4)= py.CoolProp.CoolProp.PropsSI('U','P',P(4),'T',T(4), 'Helium'); % [J/kg]
U_He(4)= U(4)*m_He_HeTank(3);

%==============================================
% Maneuver 3: critical Maneuver
%==============================================
disp("Maneuver 3");
%Helium in Propellant Tanks
%Temperature increase of Helium over Pressure regulator initial value
dT(3) = JT_coeff*(p_He_MMHTank*10^-5-P(4)*10^-5);
T_He_PR(3)=T(4)+dT(3);

for i=1:2

%MMH Tank
%Initial values
m_He_maneuver_MMHTank_final1(4)=m_He_maneuver_MMHTank(4);
m_He_maneuver_MMHTank_final2(4)=1;

while abs(m_He_maneuver_MMHTank_final1(4)-m_He_maneuver_MMHTank_final2(4))>10^-3
T_He_MMHTank(4)= (cp_He_kg*(m_He_MMHTank(3)*T_He_MMHTank(3)+m_He_maneuver_MMHTank_final1(4)*T_He_PR(3)))/(cp_He_kg*(m_He_MMHTank(3)+m_He_maneuver_MMHTank_final1(4)));
roh_He_MMHTank(4)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MMHTank,'T',T_He_MMHTank(4), 'Helium'); % [J/kg]
m_He_maneuver_MMHTank_final2(4)=V_MMH_maneuver(4)*roh_He_MMHTank(4);
if i==1
    m_He_maneuver_MMHTank_final1(4)=m_He_maneuver_MMHTank_final1(4)-0.001;
else
    m_He_maneuver_MMHTank_final1(4)=m_He_maneuver_MMHTank_final1(4)-0.001;
end

%disp(["Error", abs(m_He_maneuver_MMHTank_final1(4)-m_He_maneuver_MMHTank_final2(4))]);
if abs(m_He_maneuver_MMHTank_final1(4)-m_He_maneuver_MMHTank_final2(4))<10^-3
    break
end
end

%MON Tank
%Initial values
m_He_maneuver_MONTank_final1(4)=m_He_maneuver_MONTank(4);
m_He_maneuver_MONTank_final2(4)=1;
while abs(m_He_maneuver_MONTank_final1(4)-m_He_maneuver_MONTank_final2(4))>10^-3
T_He_MONTank(4)= (cp_He_kg*(m_He_MMHTank(3)*T_He_MONTank(3)+m_He_maneuver_MONTank_final1(4)*T_He_PR(3)))/(cp_He_kg*(m_He_MONTank(3)+m_He_maneuver_MONTank_final1(4)));
roh_He_MONTank(4)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MONTank,'T',T_He_MONTank(4), 'Helium'); % [J/kg]
m_He_maneuver_MONTank_final2(4)=V_MON_maneuver(4)*roh_He_MONTank(4);
m_He_maneuver_MONTank_final1(4)=m_He_maneuver_MONTank_final1(4)-0.001;
%disp(["Error MON Man 3", abs(m_He_maneuver_MONTank_final1(4)-m_He_maneuver_MONTank_final2(4))]);
end

%Total Helium required for maneuver 2
m_He_maneuver_final(4)= m_He_maneuver_MONTank_final2(4)+m_He_maneuver_MMHTank_final2(4);
mdot_He_maneuver_final(4)= m_He_maneuver_final(4)/t_man(4);
m_He_HeTank(4)= m_He_HeTank (3) - m_He_maneuver_final(4);
roh_He_HeTank(4) = m_He_HeTank(4)/V_tank_He ; 

%Helium in Helium Tank
U_He(5) = U_He(4)+t_man(4)*(Q_dot(3)-H(4)*mdot_He_maneuver_final(4)); % [kJ]
U(5) = U_He(5)/m_He_HeTank(4); % [kJ/kg]

P(5) = py.CoolProp.CoolProp.PropsSI('P','U',U(5),'D',roh_He_HeTank(4), 'Helium'); % [J/kg]
T(5) = py.CoolProp.CoolProp.PropsSI('T','U',U(5),'D',roh_He_HeTank(4), 'Helium'); % [J/kg]
H(5) = py.CoolProp.CoolProp.PropsSI('H','U',U(5),'D',roh_He_HeTank(4), 'Helium'); % [J/kg]


%Temperature increase of Helium over Pressure regulator
dT(3) = JT_coeff*(p_He_MMHTank*10^-5-((P(4)+P(5))/2)*10^-5);
T_He_PR(3)=(T(4)+T(5))/2+dT(3);
end

%Helium left in Propellant Tanks after third maneuver
m_He_MMHTank(4)= m_He_MMHTank(3)+m_He_maneuver_MMHTank_final2(4);
m_He_MONTank(4)=m_He_MONTank(3)+m_He_maneuver_MONTank_final2(4);

%Initial values in propellant tanks for fourth maneuver
m_He_maneuver_MMHTank(5)= V_MMH_maneuver(5)* roh_He_MMHTank(4);
m_He_maneuver_MONTank(5)=V_MON_maneuver(5)* roh_He_MONTank(4);

%==============================================
% Waiting time between Maneuver 3 and 4
%==============================================
T(6) = 293.15; % [5]
P(6) = py.CoolProp.CoolProp.PropsSI('P','T',T(6),'D',roh_He_HeTank(4), 'Helium'); % [J/kg]
H(6)= py.CoolProp.CoolProp.PropsSI('H','P',P(6),'T',T(6), 'Helium'); % [J/kg]
U(6)= py.CoolProp.CoolProp.PropsSI('U','P',P(6),'T',T(6), 'Helium'); % [J/kg]
U_He(6)= U(6)*m_He_HeTank(4);

%==============================================
% Maneuver 4
%==============================================
disp("Maneuver 4");
%Helium in Propellant Tanks
%Temperature increase of Helium over Pressure regulator initial value
dT(4) = JT_coeff*(p_He_MMHTank*10^-5-P(6)*10^-5);
T_He_PR(4)=T(6)+dT(4);

for i=1:2

%MMH Tank
%Initial values
m_He_maneuver_MMHTank_final1(5)=m_He_maneuver_MMHTank(5);
m_He_maneuver_MMHTank_final2(5)=1;

while abs(m_He_maneuver_MMHTank_final1(5)-m_He_maneuver_MMHTank_final2(5))>10^-2
T_He_MMHTank(5)= (cp_He_kg*(m_He_MMHTank(4)*T_He_MMHTank(4)+m_He_maneuver_MMHTank_final1(5)*T_He_PR(4)))/(cp_He_kg*(m_He_MMHTank(4)+m_He_maneuver_MMHTank_final1(5)));
roh_He_MMHTank(5)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MMHTank,'T',T_He_MMHTank(5), 'Helium'); % [J/kg]
m_He_maneuver_MMHTank_final2(5)=V_MMH_maneuver(5)*roh_He_MMHTank(5);
m_He_maneuver_MMHTank_final1(5)=m_He_maneuver_MMHTank_final1(5)+0.001;
end

%MON Tank
%Initial values
m_He_maneuver_MONTank_final1(5)=m_He_maneuver_MONTank(5);
m_He_maneuver_MONTank_final2(5)=1;
while abs(m_He_maneuver_MONTank_final1(5)-m_He_maneuver_MONTank_final2(5))>10^-1
T_He_MONTank(5)= (cp_He_kg*(m_He_MMHTank(4)*T_He_MONTank(4)+m_He_maneuver_MONTank_final1(5)*T_He_PR(4)))/(cp_He_kg*(m_He_MONTank(4)+m_He_maneuver_MONTank_final1(5)));
roh_He_MONTank(5)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MONTank,'T',T_He_MONTank(4), 'Helium'); % [J/kg]
m_He_maneuver_MONTank_final2(5)=V_MON_maneuver(5)*roh_He_MONTank(5);
m_He_maneuver_MONTank_final1(5)=m_He_maneuver_MONTank_final1(5)-0.001;
end

%Total Helium required for maneuver 4
m_He_maneuver_final(5)= m_He_maneuver_MONTank_final2(5)+m_He_maneuver_MMHTank_final2(5);
mdot_He_maneuver_final(5)= m_He_maneuver_final(5)/t_man(5);
m_He_HeTank (5)= m_He_HeTank (4) - m_He_maneuver_final(5);
roh_He_HeTank(5) = m_He_HeTank(5)/V_tank_He ; 

%Helium in Helium Tank
U_He(7) = U_He(6)+t_man(5)*(Q_dot(4)-H(6)*mdot_He_maneuver_final(5)); % [kJ]
U(7) = U_He(7)/m_He_HeTank(5); % [kJ/kg]

P(7) = py.CoolProp.CoolProp.PropsSI('P','U',U(7),'D',roh_He_HeTank(5), 'Helium'); % [J/kg]
T(7) = py.CoolProp.CoolProp.PropsSI('T','U',U(7),'D',roh_He_HeTank(5), 'Helium'); % [J/kg]
H(7) = py.CoolProp.CoolProp.PropsSI('H','U',U(7),'D',roh_He_HeTank(5), 'Helium'); % [J/kg]


%Temperature increase of Helium over Pressure regulator
dT(4) = JT_coeff*(p_He_MMHTank*10^-5-((P(6)+P(7))/2)*10^-5);
T_He_PR(4)=(T(6)+T(7))/2+dT(4);
end

%Helium left in Propellant Tanks after fourth maneuver
m_He_MMHTank(5)= m_He_MMHTank(4)+m_He_maneuver_MMHTank_final2(5);
m_He_MONTank(5)=m_He_MONTank(4)+m_He_maneuver_MONTank_final2(5);

%Initial values in propellant tanks for fourth maneuver
m_He_maneuver_MMHTank(6)= V_MMH_maneuver(6)* roh_He_MMHTank(5);
m_He_maneuver_MONTank(6)=V_MON_maneuver(6)* roh_He_MONTank(5);

%==============================================
% Waiting time between Maneuver 4 and 5
%==============================================
T(8) = T(7); % [5]
P(8) = py.CoolProp.CoolProp.PropsSI('P','T',T(8),'D',roh_He_HeTank(5), 'Helium'); % [J/kg]
H(8)= py.CoolProp.CoolProp.PropsSI('H','P',P(8),'T',T(8), 'Helium'); % [J/kg]
U(8)= py.CoolProp.CoolProp.PropsSI('U','P',P(8),'T',T(8), 'Helium'); % [J/kg]
U_He(8)= U(8)*m_He_HeTank(5);

%==============================================
% Maneuver 5
%==============================================
disp("Maneuver 5");
%Helium in Propellant Tanks
%Temperature increase of Helium over Pressure regulator initial value
dT(5) = JT_coeff*(p_He_MMHTank*10^-5-P(8)*10^-5);
T_He_PR(5)=T(8)+dT(5);

for i=1:2

%MMH Tank
%Initial values
m_He_maneuver_MMHTank_final1(6)=m_He_maneuver_MMHTank(6);
m_He_maneuver_MMHTank_final2(6)=1;
while abs(m_He_maneuver_MMHTank_final1(6)-m_He_maneuver_MMHTank_final2(6))>10^-1
T_He_MMHTank(6)= (cp_He_kg*(m_He_MMHTank(5)*T_He_MMHTank(5)+m_He_maneuver_MMHTank_final1(6)*T_He_PR(5)))/(cp_He_kg*(m_He_MMHTank(5)+m_He_maneuver_MMHTank_final1(6)));
roh_He_MMHTank(6)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MMHTank,'T',T_He_MMHTank(6), 'Helium'); % [J/kg]
m_He_maneuver_MMHTank_final2(6)=V_MMH_maneuver(6)*roh_He_MMHTank(6);
m_He_maneuver_MMHTank_final1(6)=m_He_maneuver_MMHTank_final1(6)+0.001;
end

%MON Tank
%Initial values
m_He_maneuver_MONTank_final1(6)=m_He_maneuver_MONTank(6);
m_He_maneuver_MONTank_final2(6)=1;
while abs(m_He_maneuver_MONTank_final1(6)-m_He_maneuver_MONTank_final2(6))>10^-1
T_He_MONTank(6)= (cp_He_kg*(m_He_MMHTank(5)*T_He_MONTank(5)+m_He_maneuver_MONTank_final1(6)*T_He_PR(5)))/(cp_He_kg*(m_He_MONTank(5)+m_He_maneuver_MONTank_final1(6)));
roh_He_MONTank(6)= py.CoolProp.CoolProp.PropsSI('D','P',p_He_MONTank,'T',T_He_MONTank(6), 'Helium'); % [J/kg]
m_He_maneuver_MONTank_final2(6)=V_MON_maneuver(6)*roh_He_MONTank(6);
m_He_maneuver_MONTank_final1(6)=m_He_maneuver_MONTank_final1(6)+0.001;
end

%Total Helium required for maneuver 5
m_He_maneuver_final(6)= m_He_maneuver_MONTank_final2(6)+m_He_maneuver_MMHTank_final2(6);
mdot_He_maneuver_final(6)= m_He_maneuver_final(6)/t_man(6);
m_He_HeTank (6)= m_He_HeTank (5) - m_He_maneuver_final(6);
roh_He_HeTank(6) = m_He_HeTank(6)/V_tank_He ; 

%Helium in Helium Tank
U_He(9) = U_He(8)+t_man(6)*(Q_dot(5)-H(8)*mdot_He_maneuver_final(6)); % [kJ]
U(9) = U_He(9)/m_He_HeTank(6); % [kJ/kg]

P(9) = py.CoolProp.CoolProp.PropsSI('P','U',U(9),'D',roh_He_HeTank(6), 'Helium'); % [J/kg]
T(9) = py.CoolProp.CoolProp.PropsSI('T','U',U(9),'D',roh_He_HeTank(6), 'Helium'); % [J/kg]
H(9) = py.CoolProp.CoolProp.PropsSI('H','U',U(9),'D',roh_He_HeTank(6), 'Helium'); % [J/kg]


%Temperature increase of Helium over Pressure regulator
dT(5) = JT_coeff*(p_He_MMHTank*10^-5-((P(8)+P(9))/2)*10^-5);
T_He_PR(5)=(T(8)+T(9))/2+dT(5);
end

%Helium left in Propellant Tanks after fifth maneuver
m_He_MMHTank(6)= m_He_MMHTank(5)+m_He_maneuver_MMHTank_final2(6);
m_He_MONTank(6)=m_He_MONTank(5)+m_He_maneuver_MONTank_final2(6);

%==============================================
% Results
%==============================================
% Total Mass of Helium required
m_He_total= m_He_maneuver_final(1)+m_He_maneuver_final(2)+m_He_maneuver_final(3)+m_He_maneuver_final(4)+m_He_maneuver_final(5)+m_He_maneuver_final(6);
mdot_He_maneuver_final= m_He_maneuver_final/t_man;
%Total Helium Mass at the beginning of operation in Propellant Tanks and
%Helium Tank
disp(["Helium per maneuver MON Tank",m_He_maneuver_MONTank_final2])
disp(["Helium per maneuver MMH Tank",m_He_maneuver_MMHTank_final2])


m_He_allTanks= m_He_HeTank(1)+m_He_MMHTank(1)+m_He_MONTank(1);

disp("Total Helium");
disp(["Total Helium Mass in all tanks[kg]",m_He_allTanks]);

disp("Helium Tank");
disp(["Helium Mass[kg]",m_He_HeTank]);
disp(["Total Helium Mass Maneuvers [kg]",m_He_total]);

disp("Helium MMH Tank");
disp(["Ullage He in MMH[kg]",m_He_MMHTank(1)])
disp(["Helium Mass in MMH tank[kg]",m_He_MMHTank]);

disp("Helium in MON Tank");
disp(["Ullage[kg]",m_He_MONTank(1)])
disp(["Helium Mass[kg]",m_He_MONTank]);

disp(["Helium mass flow each maneuver[kg]",mdot_He_maneuver_final]);

disp(["He Pressure[bar]", P*(1e-5)]);
disp(["He Temperature[ºC]", T-273.15]);
disp(["He Density [kg/m^3]",roh_He_HeTank]);
disp(["Mass of Helium required each maneuver[kg]",m_He_maneuver_final]);
disp(["Start Iteration",m_He_maneuver_MMHTank_final1]);
disp(["Finished Iteration",m_He_maneuver_MMHTank_final2]);
disp(["Mass of Helium required each maneuver MON Tank[kg]",m_He_maneuver_MONTank_final2]);

fig1=figure(1);
plot(T-273.15);
xlabel('Point in time');
ylabel('Temperature [ºC]');
title('Temperature');

fig2=figure(2);
plot(P*(1e-5));
xlabel('Point in time');
ylabel('Pressure [bar]');
title('Pressure');

fig3=figure(3);
plot(roh_He_HeTank);
xlabel('Point in time');
ylabel('Density [kg/m^3]');
title('Density');

fig4=figure(4);
plot(t_man_tot/60,P*(1e-5));
xlabel ('Total Maneuver time[min]');
ylabel('Pressure [bar]');

fig5=figure(5);
plot(t_man_tot/60,T-273.15);
xlabel ('Total Maneuver time[min]');
ylabel('Temperature [ºC]');

fig6=figure(6);
plot(t_man_tot0/60,roh_He_HeTank);
xlabel ('Total Maneuver time[min]');
ylabel('Density [kg/m^3]');

fig8=figure(8);
plot(m_He_HeTank);
xlabel ('Point in time');
ylabel('Mass Helium in Helium Tank[kg]');

fig9=figure(9);
plot(m_He_MMHTank);
xlabel ('Point in time');
ylabel('Mass Helium in MMH Tank[kg]');

fig10=figure(10);
plot(m_He_MONTank);
xlabel ('Point in time');
ylabel('Mass Helium in MON Tank[kg]');