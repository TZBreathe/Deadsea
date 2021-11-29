
%
clear;

load BrOcv;

SoC_LUT = BrOcv.Dims.soc';   % SoC
OCV_LUT = BrOcv.Components.ocv(:,5); % OCV, V


R0 = 0.0005; 
R1 = 0.0020; 
C1 = 250000; 
Q0=10;

a = R1/15; % Power fade coefficient, A^-0.5h^-0.5
Tamb = 25; % Ambient temperature, deg
I_ch = 10; % Charge current, V
I_disch = -10; % Discharge current, V
V_max = 4.2; % upper voltage cut-off, V
V_min = 3; % lower voltage cut-off, V
t_step = 1; % Time step,s 
I = 10; % Charge/Discharge Current Magnitude, A
%% Initial Conditions
SoC_ini = 1; % Initial state-of-charge
Em_ini = interp1(SoC_LUT,OCV_LUT,SoC_ini); % Intial OCV
Vcell_ini = Em_ini; % Initial cell voltage
Q=Q0; % Actual cell capacity at beginging
i = 1;
Em(1) = Em_ini;
Vcell(1) = Vcell_ini;
SoC(1) = SoC_ini; % SoC as defined by nominal capacity
SoC_scaled(i)=SoC_ini; % Scaled SoC as defined by actual cell capacity
t(1)=0;
%% Simulation
for k = 1:500
 
 
%  while Vcell(i) < V_max
%  % Update time
%  t = t + t_step;
%  i = i + 1;
%  
%  % Update battery state
% % SoC as defined by nominal capacity 
%  SoC(i)=SoC(i-1)+I_ch*t_step/3600/Q0;
% 
% % Look up OCV using scaled SoC
%  Em(i)= interp1(SoC_LUT,OCV_LUT,SoC(i)); 
%  % Calculate cell voltage
%  Vcell(i)=Em(i)+I_ch*(R1*(1-exp(-t/(R1*C1)))+R0);
%  end
%  

 time_start(k)=t;
 SoC(i)=SoC_ini;
 Em(i)= interp1(SoC_LUT,OCV_LUT,SoC(i)); 
 Vcell(i)=Em(i);
 while Vcell(i) > V_min
 
     
 % Update time 
 t= t + t_step;
 i = i + 1;

 
 % Upddate battery state 
 SoC(i) = SoC(i-1)+I_disch*t_step/3600/Q0;
 Em(i) = interp1(SoC_LUT,OCV_LUT,SoC(i));
 % Calculate Cell voltage
 Vcell(i) = Em(i) + I_disch*(R1*(1 - exp(-t/(R1*C1))) + R0);
 end
 
 time_end(k)=t;
 
 % Update cell capacity afte each charge/discharge
 R1 = R1 + a;
 
 
end
%% Plot
CapList=(time_end-time_start)*I_ch/3600;
figure()
plot(CapList);
xlabel('Cycle');ylabel('Capacity/Ah');
figure();
plot(linspace(0,time_end(2)-time_start(2),time_end(2)-time_start(2)+1),Vcell(time_start(2):time_end(2)));
hold on
plot(linspace(0,time_end(199)-time_start(199),time_end(199)-time_start(199)+1),Vcell(time_start(199):time_end(199)));
hold off
legend('2nd cycle','500th cycle')
xlabel('Time');ylabel('Voltage');
