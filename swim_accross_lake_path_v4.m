clc; clear all; close all;

%% Definitions
h = 1;
u = 3;
x = linspace(0,h,1000);
delta_x = x(2)-x(1);
for i = 1:length(x)
    v(i) = wave_speed(x(i));
end
v_max = max(v);
v_min = min(v);

%% Path calculation
lambda = linspace(-0.99/(u+v_max),0.99/(u-v_min),9); %-1/(u+v_max) < lambda < 1/(u-v_min)
for j = 1:length(lambda)
    y(1,j) = 0; %IC
    T(1,j) = 0; %IC
    for i = 1:length(x)
        y(i+1,j) = y(i,j) + dydx(u,lambda(j),x(i))*delta_x; %y path
        T(i+1,j) = T(i,j) + dtdx(u,lambda(j),x(i))*delta_x; %time taken along path
    end
end
y(end,:) = [];
T(end,:) = [];

%% Plots
figure(1) %wave speed
plot(x,v);
title('water speed upwards')
xlabel('river width (km)')
ylabel('water speed (km/hr)')
set(gcf,'color','w');

figure(2) %path
for j = 1:length(lambda)
    if j == 7
        plot(x,y(:,7),'DisplayName',strcat('\lambda = ',num2str(lambda(7)),' and T = ',num2str(T(end,7))),'LineWidth',2,'color',[0.9290 0.6940 0.1250])
    else
        plot(x,y(:,j),'DisplayName',strcat('\lambda = ',num2str(lambda(j)),' and T = ',num2str(T(end,j))),'LineWidth',2)
    end
    hold on
end
legend
grid on
title('Optimal path across river')
xlabel('river width (km)')
ylabel('river length (km)')
set(gcf,'color','w');

figure(3)
plot(x,y(:,7),'DisplayName',strcat('\lambda = ',num2str(lambda(7)),' and T = ',num2str(T(end,7))),'LineWidth',2,'color',[0.9290 0.6940 0.1250])
grid on
% title('Optimal path across river')
title(strcat('Optimal path across river: \lambda = ',num2str(lambda(7))))
xlabel('river width (km)')
ylabel('river length (km)')
ylim([-1,2.5])
set(gcf,'color','w');
%% Functions
function f = dydx(u,lambda,x)
    %water speed function
    v = wave_speed(x);
    
    %integral equation;
    S = - u*lambda / (1+lambda*v);
    C = cos(asin(S));
    f = (u*S+v)/(u*C); %dy/dx = f = Vy/Vx
end

function f = dtdx(u,lambda,x)
    v = wave_speed(x);
    C = sqrt((1+lambda*v)^2-(u*lambda)^2)/(1+lambda*v);
    
    f = 1/(u*C);
end

function v = wave_speed(x)
    v = 3*x;
end