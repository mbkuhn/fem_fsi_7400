function [ksi_G,eta_G,W_G] = tri_Gpts_3()
% Gauss points for linear triangle element, ng = 3
% No inputs
% Outputs: ksi_G coordinates, eta_G coordinates, W_G weights

ksi_G = [1/6 1/6 2/3];
eta_G = [1/6 2/3 1/6];
W_G   = 1/6*ones(1,3);
end