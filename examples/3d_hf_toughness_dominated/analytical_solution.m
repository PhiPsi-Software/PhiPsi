% Written By: Shi Fang, 2019-11-06
% Website: phipsi.top

% Analysis solution for the 3D penny-shaped fracturing model.
clear all; close all; clc; format compact;  format long;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Input parameters.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Gupta_2018_Coupled hydromechanical-fracture simulations of nonplanar three-dimensional hydraulic fracture propagation_Section 6.2
% E            = 17e9;        
% v            = 0.2;         % Poisson's ratio
% K_IC         = 1.46e6;
% u            = 0.1;         % viscosity
% u            = 0.0001;      % viscosity
% Q_0          = 0.001;       % Rate of flow.
% Total_time   = 15;          % fracturing time(s)

% Salimzadeh_2017_Finite element simulations of interactions between multiple hydraulic fractures in a poroelastic rock_Section 5.2
E            = 17e9;        
v            = 0.25;       % Poisson's ratio
K_IC         = 2.0e6;
% u            = 0.1;         % viscosity
u            = 1.0e-4;     % viscosity
Q_0          = 0.01;       % Rate of flow.
Total_time   = 20;         % fracturing time(s)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot control parameters.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot control(1: plot R vs time),(2: plot injection pressure vs time),(3: plot w vs R),(4: plot pressure vs R)
key_Plot(1:4)       = [1,  1,  1,  1];
Number_Plot_Points  = 500;    % number of plot points

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Main program.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Intermediate variable.
E_1 = E/(1-v^2);
K_1 = 4*sqrt(2/pi)*K_IC;
u_1 = 12*u;
disp(['Visocity is    ', num2str(u ,'%0.4f'),' Pa.s.'])
disp(['Total time is  ', num2str(Total_time ,'%0.4f'),' s.'])

%Check fracturing type.
if u>0.01
    Fracturing_Type=1; 
	disp(['*************************************************'])
    disp(['   Warnning :: viscosity-dominated fracturing'])
	disp(['*************************************************'])
elseif u<=0.01
    Fracturing_Type=2; 
	disp(['*************************************************'])
    disp(['Warnning :: toughness-dominated fracturing'])
	disp(['*************************************************'])
end

%####################################
%  viscosity-dominated fracturing
%####################################
if Fracturing_Type==1 
	%******************************
	%       plot R vs time
	%******************************
    if key_Plot(1)==1
	    disp(['Plotting R vs time curve......'])
	    for i =1:Number_Plot_Points+1
		    t(i) = 0.0 + (i-1)*Total_time/Number_Plot_Points;
			R(i) = 0.6955*(Q_0^3*E_1/u_1)^(1.0/9.0)*(t(i)^(4/9));
		end
		figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on
		plot(t,R)
		title('\it R vs time','FontName','Times New Roman','FontSize',16)
		hold on		
	end
	% Save data
	disp(['Saving data......'])
	Text_t = t';
	Text_R = R';	 
	save Viscosity-Dominated_time_values.txt Text_t -ascii
    save Viscosity-Dominated_R_values.txt Text_R -ascii
	
	%************************************
	%  plot injection pressure vs time
	%************************************
    if key_Plot(2)==1
	    disp(['Plotting injection pressure vs time curve......'])
	    for i =1:Number_Plot_Points+1
		    t(i) = 0.0 + (i-1)*Total_time/Number_Plot_Points;
			R(i) = 0.6955*(Q_0^3*E_1/u_1)^(1.0/9.0)*(t(i)^(4/9));
			injection_P(i) = (0.8593-0.2387-0.09269*(log(0.01/R(i))))*(u_1*E_1^2)^(1.0/3.0)*(t(i))^(-1.0/3.0);
		end
		figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on
		plot(t,injection_P)
		title('\it injection pressure vs time','FontName','Times New Roman','FontSize',16)
		hold on		
	end
	% Save data
	disp(['Saving data......'])
	Text_injection_P = injection_P';	 
    save Viscosity-Dominated_injection_Pressure_values.txt Text_injection_P -ascii
	
	%************************************
	%  plot the final aperture (w vs R)
	%************************************
    if key_Plot(3)==1
	    disp(['Plotting the final aperture (w vs R)......'])
		R_final = 0.6955*(Q_0^3*E_1/u_1)^(1.0/9.0)*(Total_time^(4/9));
	    for i =1:Number_Plot_Points+1
			Final_r_points(i) = 0.0 +  (i-1)*R_final/Number_Plot_Points;
			tem   = Final_r_points(i)/R_final;
			
			part1 = 1.9093*(1.0-tem)^(2.0/3.0);
			part2 = 0.0705*(13.0*tem - 6.0)*((1.0-tem)^(2.0/3.0));
			
			part3 = 0.236*(1.0-tem^2)^0.5;
			part4 = 0.236*tem*acos(tem);
			part5 = ((Q_0^3*u_1^2)/(E_1^2))^(1.0/9.0)*Total_time^(1.0/9.0);
			w(i) = 0.6955*(part1 + part2 + part3 - part4)*part5*1000.0;
		end
		
		figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on
		plot(Final_r_points,w)
		title('\it aperture(mm) vs R','FontName','Times New Roman','FontSize',16)
		hold on		
	end
	% Save data
	disp(['Saving data......'])
	Text_Final_r_points = Final_r_points';	 
    save Viscosity-Dominated_Final_r_points_values.txt Text_Final_r_points -ascii	
	Text_w = w';	 
    save Viscosity-Dominated_Final_w.txt Text_w -ascii		
	
	%************************************
	%  plot the final pressure (p vs R)
	%************************************
    if key_Plot(4)==1
	    disp(['Plotting the final pressure (p vs R)......'])
		R_final = 0.6955*(Q_0^3*E_1/u_1)^(1.0/9.0)*(Total_time^(4/9));
	    for i =1:Number_Plot_Points+1
			Final_r_points(i) = 0.0 +  (i-1)*R_final/Number_Plot_Points;
			tem   = Final_r_points(i)/R_final;
			p(i) = (0.8593-0.2387/((1.0-tem)^(1/3))-0.09269*log(tem))*(u_1*E_1^2)^(1.0/3.0)*(Total_time)^(-1.0/3.0)/1.0e6;
		end
		
		
		figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on
		plot(Final_r_points,p)
		title('\it pressure(MPa) vs R','FontName','Times New Roman','FontSize',16)
		hold on		
	end
	% Save data
	disp(['Saving data......'])
	Text_p = p';	 
    save Viscosity-Dominated_Final_p.txt Text_p -ascii		
	
%####################################
%   toughness-dominated fracturing
%####################################
elseif Fracturing_Type==2
%******************************
	%       plot R vs time
	%******************************
    if key_Plot(1)==1
	    disp(['Plotting R vs time curve......'])
	    for i =1:Number_Plot_Points+1
		    t(i) = 0.0 + (i-1)*Total_time/Number_Plot_Points;
			Mu   = u_1*((Q_0^3*E_1^13)/(K_1^18* t(i)^2))^(1.0/5.0);  %Savitski_2002_Eq.(34)
			R(i) = (0.8546+Mu*(-0.7349))*((Q_0^2*E_1^2*t(i)^2)/K_1^2)^(1.0/5.0);
		end
		figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on
		plot(t,R)
		title('\it R vs time','FontName','Times New Roman','FontSize',16)
		hold on		
	end
	% Save data
	disp(['Saving data......'])
	Text_t = t';
	Text_R = R';	 
	save toughness-Dominated_time_values.txt Text_t -ascii
    save toughness-Dominated_R_values.txt Text_R -ascii
	
	%************************************
	%  plot injection pressure vs time
	%************************************
    if key_Plot(2)==1
	    disp(['Plotting injection pressure vs time curve......'])
	    for i =1:Number_Plot_Points+1
		    t(i)           =  0.0 + (i-1)*Total_time/Number_Plot_Points;
			Mu             =  u_1*((Q_0^3*E_1^13)/(K_1^18*t(i)^2))^(1.0/5.0);  %Savitski_2002_Eq.(34)
			Epsilon        = (K_1^6/(E_1^6*Q_0*t(i)))^(1.0/5.0);                %Savitski_2002_Eq.(32-1)
			R(i)           = (0.8546+Mu*(-0.7349))*((Q_0^2*E_1^2*t(i)^2)/K_1^2)^(1.0/5.0);
			rho            =  0.01/R(i);                                        % ρ=r/R, r is selected as 0.01 meter.
			PI             = (0.3004+Mu*(0.6380-1.709*(1/3*log(rho)-1/5*log(1.0-rho^2)))); %Savitski_2002_Eqs.(74),(84),(88)
			injection_P(i) =  E_1*Epsilon*PI;                                   %Savitski_2002_Eq.(13)
		end
		figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on
		plot(t,injection_P)
		title('\it injection pressure vs time','FontName','Times New Roman','FontSize',16)
		hold on		
	end
	% Save data
	disp(['Saving data......'])
	Text_injection_P = injection_P';	 
    save toughness-Dominated_injection_Pressure_values.txt Text_injection_P -ascii
	
	
	%************************************
	%  plot the final aperture (w vs R)
	%************************************
    if key_Plot(3)==1
	    disp(['Plotting the final aperture (w vs R)......'])
		Mu      =  u_1*((Q_0^3*E_1^13)/(K_1^18*Total_time^2))^(1.0/5.0);  %Savitski_2002_Eq.(34)		
		R_final = (0.8546+Mu*(-0.7349))*((Q_0^2*E_1^2*Total_time^2)/K_1^2)^(1.0/5.0);
		Epsilon = (K_1^6/(E_1^6*Q_0*Total_time))^(1.0/5.0);               %Savitski_2002_Eq.(32-1)		
		L       = ((Q_0^2*E_1^2*t(i)^2)/K_1^2)^(1.0/5.0);                     %Savitski_2002_Eq.(32-2)		
	    for i =1:Number_Plot_Points+1
			Final_r_points(i) = 0.0 +  (i-1)*R_final/Number_Plot_Points;
			tem   = Final_r_points(i)/R_final;
			omega_k0 = (3/8/pi)^(1/5)*(1-tem^2)^(1/2);
			%numerical integration
            integration_sum   = 0;	
			integration_point = 100; 
			for j= 1:integration_point
			   delta_x =  (1.0-tem)/integration_point;
			   x =  tem + j*delta_x;
			   integration_sum = integration_sum + delta_x*sqrt((1-x^2)/(x^2-tem^2))*asin(x);
			end
			% integration_sum
			
			omega_k1 = 0.8264-1.23972*(-0.10685*(1-tem^2)^0.5 + tem*acos(tem)-6/5*integration_sum);
			w(i) = Epsilon*L*(omega_k0 + Mu*omega_k1)*1000; %Savitski_2002_Eq.(12)	
		end
		figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on
		plot(Final_r_points,w)
		title('\it aperture(mm) vs R','FontName','Times New Roman','FontSize',16)
		hold on		
	end
	% Save data
	disp(['Saving data......'])
	Text_Final_r_points = Final_r_points';	 
    save toughness-Dominated_Final_r_points_values.txt Text_Final_r_points -ascii	
	Text_w = w';	 
    save toughness-Dominated_Final_w.txt Text_w -ascii		
	
	%************************************
	%  plot the final pressure (p vs R)
	%************************************
    if key_Plot(4)==1
	    disp(['Plotting the final pressure (p vs R)......'])
		Mu      =  u_1*((Q_0^3*E_1^13)/(K_1^18*Total_time^2))^(1.0/5.0);  %Savitski_2002_Eq.(34)		
		R_final = (0.8546+Mu*(-0.7349))*((Q_0^2*E_1^2*Total_time^2)/K_1^2)^(1.0/5.0);
		Epsilon = (K_1^6/(E_1^6*Q_0*Total_time))^(1.0/5.0);               %Savitski_2002_Eq.(32-1)		
		L       = ((Q_0^2*E_1^2*t(i)^2)/K_1^2)^(1.0/5.0);                     %Savitski_2002_Eq.(32-2)	
		
	    for i =1:Number_Plot_Points+1
			Final_r_points(i) = 0.0 +  (i-1)*R_final/Number_Plot_Points;
			tem   = Final_r_points(i)/R_final;
			PI    = (0.3004+Mu*(0.638-1.709*(1/3*log(tem)-1/5*log(1.0-tem^2)))); %Savitski_2002_Eqs.(74),(84),(88)
			p(i) =  E_1*Epsilon*PI/1.0e6;                                   %Savitski_2002_Eq.(13)
		end
		
	
	
		figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on
		plot(Final_r_points,p)
		title('\it pressure(MPa) vs R','FontName','Times New Roman','FontSize',16)
		hold on		
	end
	% Save data
	disp(['Saving data......'])
	Text_p = p';	 
    save toughness-Dominated_Final_p.txt Text_p -ascii	
end
disp(['All done!'])