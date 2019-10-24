clear, clc, close all

% load parameter set structure
load('base_param_set')

% set max slider values, and parameter units
[sliderMax, paramUnit] = gui_config;

% instantiate gui object
gui = CalciumGUI(P, 'gui', sliderMax, paramUnit);







