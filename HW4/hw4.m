%%
clc; clear all; close all;

t0 = 0;
data = importdata('homerun_data_HW4.txt');
tspan = data(:, 1);
rho = data(:, 2);
alpha = data(:, 3)*pi/180;
beta = data(:, 4)*pi/180;

i = 1;
t_i_1 = t0;
Xstar_i_1 = data(1, :);
xbar_0 = [0; 0; 0]; % COMBAK: temporarily setting the apriori state data to 0's but I dont think that's correct
xhat_i_1 = xbar_0;
P_i_1 = Pbar_0;
