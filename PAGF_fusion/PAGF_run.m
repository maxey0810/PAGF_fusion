close all;clc;clear;
addpath metrics\

% I\O setting
input_path = 'input/';

img_name = 'Camp1827.jpg'; 
IR = im2double(imread(strcat(input_path,'IR/', img_name)));

VI = im2double(imread(strcat(input_path,'VI/', img_name)));

% If the energy information in visible images becomes unreliable under extreme conditions, 
% it is recommended to manually set a large nonlinear weight parameter (e.g., >20 in 'meting.jpg').

   F=PAGF_fusion(IR,VI); %adaptive parameters
   %F=PAGF_fusion(IR,VI,b_fixed);%manually set a large nonlinear weight parameter

figure;imshow([IR,VI,F]);title('PAGF');

% metrics
F=im2uint8(F);
IR=im2uint8(IR);
VI=im2uint8(VI);

PAGF_CE = metricsCross_entropy(IR,VI,F);
PAGF_AG = metricsAvg_gradient(IR,VI,F);
PAGF_SF = metricsSpatial_frequency(IR,VI,F);
PAGF_SCD=metricsScd(IR,VI,F);
PAGF_EI = metricsEdge_intensity(IR,VI,F);
PAGF_VIF = metricsVif(IR,VI,F);




