clc; clear all; close all;

%% Calculations

%create empty history file
history= struct('y',[],'fval',[]); % create the structure with empty fields
history(1)= []; % create the empty structure array with the specified fields
save('history.mat');

%constants
h = 1; %width of river
u = 3; %speed of swimmer (by himself)
n = 40; %number of points
x = linspace(0,h,n); %horizontal distance across river

%enter in start and end positions
y_start = 0;
y_end = 0;

objective = @(y) objective_fun(u,h,x,y);
x0 = [y_start zeros(1,n-2) y_end];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [y_start -Inf(1,n-2) y_end];
ub = [y_start Inf(1,n-2) y_end];
nonlcon = [];
options = optimoptions(@fmincon,'OutputFcn',@outfun2);

[y_optimized,fval,~,output] = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

set(gcf, 'Position',  [100, 100, 1600, 900], 'color', 'w')
plot(x,y_optimized,'-o') %plot of optimal path
title('Optimal path using solver directly')


%% Video

load('history.mat');
[m, num_iterations] = size(history);
figure()
set(gcf,'renderer','painters')
for i = 1:num_iterations
y_vals = history(i).y; %each iteration of y

    set(gcf, 'Position',  [100, 100, 1600, 900], 'color', 'w')
    plot(x,y_vals,'-o')
    grid on
    xlabel('river width (km)')
    ylabel('river height (km)')
    xlim([0 h])
    ylim([-0.3 0.3]);

    movieVector(i) = getframe(gcf);
end

myWriter = VideoWriter('Math_puzzle_swim_across_lake_v_3x');
myWriter.FrameRate = 5;
open(myWriter);
writeVideo(myWriter, movieVector)
close(myWriter)

%% Functions

%For movie: saves iterates in separate history file
function stop = outfun2(y,optimValues,state)
    stop = false;
    load('history.mat','history');
    history(end+1)= struct('y',y,'fval',optimValues.fval);
    save('history.mat');
end

%Defining combined error between points and curve
function T = objective_fun(u,h,x,y)

    %Define current velocity
    for i = 1:length(x)
        v(i) = 3*x(i);
    end

    %Find time between points of curve
    dx = x(2) - x(1);
    for i = 1:length(x)-1
        dy = y(i+1)-y(i);
        %The following expressions is derived using these 3 discrete formulas
        % u_x*dt = dx; %horizontal velocity
        % u_y*dt + v*dt = dy; %vertical velocity
        % u^2 = u_x^2 + u_y^2; %pythagorus
%         dt(i) = (-2*v(i)*dy + sqrt((2*v(i)*dy)^2+4*(u^2-v(i)^2)*(dy^2+dx^2)))/(2*u^2-2*v(i)^2);
        dt(i) = (-v(i)*dy + sqrt((v(i)*dy)^2+(u^2-v(i)^2)*(dy^2+dx^2)))/(u^2-v(i)^2);
    end

    T = sum(dt); %total time taken
end