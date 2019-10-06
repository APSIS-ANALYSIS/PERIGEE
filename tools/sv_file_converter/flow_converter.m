% data will load the .flow file and fourier fit in matlab will generate
% the Fourier series for the inflow data. We use 6 mode by defaut. but 
% make sure the fit is not over interpolated.
clear all; clc;

data = load('/home/jliu/SV-Examples/PH_control_2017_400_PCMRI_rigid/2017_400_PCMRI.flow');

f = fit(data(:,1), -data(:,2), 'fourier7');

plot(f, data(:,1), -data(:,2));

num_mode = 7;

file = '/home/jliu/inflow_fourier_series.txt';

fid = fopen(file, 'w');

fprintf(fid, '# num_of_mode w perid\n');
fprintf(fid, '%d %e %e\n', num_mode, f.w, data(end,1) );

fprintf(fid, '\n');
fprintf(fid, '# coef_a with length num_of_mode + 1\n');

fprintf(fid, '%e\t', f.a0);
fprintf(fid, '%e\t', f.a1);
fprintf(fid, '%e\t', f.a2);
fprintf(fid, '%e\t', f.a3);
fprintf(fid, '%e\t', f.a4);
fprintf(fid, '%e\t', f.a5);
fprintf(fid, '%e\t', f.a6);
fprintf(fid, '%e\n', f.a7);

fprintf(fid, '\n');
fprintf(fid, '# coef_b with length num_of_mode + 1\n');

fprintf(fid, '%e\t', 0.0);
fprintf(fid, '%e\t', f.b1);
fprintf(fid, '%e\t', f.b2);
fprintf(fid, '%e\t', f.b3);
fprintf(fid, '%e\t', f.b4);
fprintf(fid, '%e\t', f.b5);
fprintf(fid, '%e\t', f.b6);
fprintf(fid, '%e\n', f.b7);

fclose(fid);

% EOF