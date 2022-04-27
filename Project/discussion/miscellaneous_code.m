%%
clc; clear all; close all;

%% Finding inv_M for a Hex
syms Ct l ly lx Cq real
Ct = 8.54858e-06;
Cq = 1.2246e-07;
L = 0.5477666432;
Lx = 0.4743798283;
Ly = 0.2738833216;
    
M = [Ct Ct Ct Ct Ct Ct;
     -Ct*L Ct*L Ct*Ly -Ct*Ly -Ct*Ly Ct*Ly;
     0 0 -Ct*Lx Ct*Lx -Ct*Lx Ct*Lx;
     Cq -Cq Cq -Cq -Cq Cq]

inv_M = pinv(M)

%%
Ct = 5.84e-6;
Cq = 1.4521e-6;
Lrf = 0.22;
Lrb = 0.2;
Lp = 0.13;

inv_M = [1/(4*Ct), -1/(2*Ct*(Lrb + Lrf)), -1/(4*Ct*Lp), -Lrb/(2*Cq*(Lrb + Lrf));
         1/(4*Ct),  1/(2*Ct*(Lrb + Lrf)),  1/(4*Ct*Lp), -Lrf/(2*Cq*(Lrb + Lrf));
         1/(4*Ct),  1/(2*Ct*(Lrb + Lrf)), -1/(4*Ct*Lp),  Lrb/(2*Cq*(Lrb + Lrf));
         1/(4*Ct), -1/(2*Ct*(Lrb + Lrf)),  1/(4*Ct*Lp),  Lrf/(2*Cq*(Lrb + Lrf))]