% Written By: Shi Fang, 2017
% Website: phipsi.top
% Email: shifang@hyit.edu.cn

% Analysis solution for the plane strain KGD fracturing model.

clear all; close all; clc; format compact;  format long;

E    = 20e9;        
v    = 0.2;         
K_IC = 0.1e6;
u    = 0.1;         
Q_0  = 0.001;      
time = 75.00;

% Intermediate variable.
E_1 = E/(1-v^2);
K_1 = 4*sqrt(2/pi)*K_IC;
u_1 = 12*u;
k = K_1/(((E_1^3)*u_1*Q_0)^0.25);
km = 4*sqrt((2/pi))*(K_IC*(1-v^2)/E)*(E/(12*u*Q_0*(1-v^2)))^0.25;

disp(['k is ', num2str(k ,'%0.4f'),'.'])

% Calculate data.
for i=1:time*10
    t(i)=(i-1)*(1/10);
	
	L(i)=(E_1*(Q_0^3)*(t(i)^4)/u_1)^(1/6); 
	
    l(i)=0.616*L(i);
	W_max(i) = ((u_1/E_1/t(i))^(1/3))*0.616*L(i)*1.84;

	P_max(i) =0.547* E/(1-v^2)*(12*u*(1-v^2)/E/t(i))^(1/3);
end

% New figure.
figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
plot(t,l);
title('\it t vs L','FontName','Times New Roman','FontSize',16)

% New figure.
figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
hold on
plot(t,W_max)
title('\it t vs W','FontName','Times New Roman','FontSize',16)

% New figure.
figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
plot(t,P_max)
title('\it t vs borehole pressure','FontName','Times New Roman','FontSize',16)
