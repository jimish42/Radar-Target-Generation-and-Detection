% TODO : Find the Bsweep of chirp for 1 m resolution
c = 3*10^8;
delta_r = 1;
b_sweep = c / 2 * delta_r;

% TODO : Calculate the chirp time based on the Radar's Max Range
r_max = 300;
% 5.5 time of the trip time for max range
ts = 5.5 * (r_max * 2 / c);

% TODO : define the frequency shifts 
beat_f = [0 1.1e6 13e6 24e6];
calc_r = c * ts * beat_f / (2 * b_sweep);

% Display the calculated range
disp(calc_r);
