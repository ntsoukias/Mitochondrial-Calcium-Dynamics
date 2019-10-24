function [sliderMax, paramUnit] = gui_config

% Configuration for GUI
% Specify slider maximum values, and parameter units

k = 1;
% volCt 
sliderMax(k, 1) = 1; 
paramUnit{k, 1} = 'pL';
k = k + 1;
% volER 
sliderMax(k, 1) = 1; 
paramUnit{k, 1} = 'pL';
k = k + 1;
% volMt
sliderMax(k, 1) = 1; 
paramUnit{k, 1} = 'pL';
k = k + 1;
% Vip3r 
sliderMax(k, 1) = 4; 
paramUnit{k, 1} = '1/s';
k = k + 1;
% Vserca 
sliderMax(k, 1) = 100; 
paramUnit{k, 1} = 'uM/s';
k = k + 1;
% kserca 
sliderMax(k, 1) = 1; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% ip3 
sliderMax(k, 1) = 4; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% Vmcu 
sliderMax(k, 1) = 15; 
paramUnit{k, 1} = 'uM/s';
k = k + 1;
% kmcu 
sliderMax(k, 1) = 10; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% Vncx 
sliderMax(k, 1) = 150; 
paramUnit{k, 1} = 'uM/s';
k = k + 1;
% kncx 
sliderMax(k, 1) = 80; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% kna 
sliderMax(k, 1) = 15; 
paramUnit{k, 1} = 'mM';
k = k + 1;
% N 
sliderMax(k, 1) = 15; 
paramUnit{k, 1} = 'mM';
k = k + 1;
% N_u 
sliderMax(k, 1) = 15; 
paramUnit{k, 1} = 'mM';
k = k + 1;
% leak_e_u 
sliderMax(k, 1) = 0.1; 
paramUnit{k, 1} = '1/s';
k = k + 1;
% leak_e_c 
sliderMax(k, 1) = 0.1; 
paramUnit{k, 1} = '1/s';
k = k + 1;
% leak_u_c 
sliderMax(k, 1) = 0.1; 
paramUnit{k, 1} = '1/s';
k = k + 1;
% cI
sliderMax(k, 1) = 1; 
paramUnit{k, 1} = '';
k = k + 1;
% cS
sliderMax(k, 1) = 1; 
paramUnit{k, 1} = '';
k = k + 1;
% cM
sliderMax(k, 1) = 1; 
paramUnit{k, 1} = '';
k = k + 1;
% cN
sliderMax(k, 1) = 1; 
paramUnit{k, 1} = '';
k = k + 1;
% bt_c 
sliderMax(k, 1) = 800; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% K_c
sliderMax(k, 1) = 20; 
paramUnit{k, 1} = '';
k = k + 1;
% bt_e 
sliderMax(k, 1) = 500000; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% K_e
sliderMax(k, 1) = 1000; 
paramUnit{k, 1} = '';
k = k + 1;
% bt_m 
sliderMax(k, 1) = 1000000; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% K_m
sliderMax(k, 1) = 800; 
paramUnit{k, 1} = '';
k = k + 1;
% bt_u 
sliderMax(k, 1) = 800; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% K_u
sliderMax(k, 1) = 30; 
paramUnit{k, 1} = '';
k = k + 1;
% a2 
sliderMax(k, 1) = 0.5; 
paramUnit{k, 1} = '1/(uM*s)';
k = k + 1;
% d1 
sliderMax(k, 1) = 0.3; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% d2 
sliderMax(k, 1) = 3; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% d3 
sliderMax(k, 1) = 4; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% d5 
sliderMax(k, 1) = 0.6; 
paramUnit{k, 1} = 'uM';
k = k + 1;
% distance
sliderMax(k, 1) = 300; 
paramUnit{k, 1} = 'nm';

