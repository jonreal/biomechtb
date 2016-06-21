M = 50;     % operating torque
Un = 18;    % nominal motor voltage
n = 15000;
n0 = 31000; % no load speed @ Un
delta_n_div_delta_M = 48.4;  % speed/torque gradient of the motor

Vcc = (Un/n0 * (n + delta_n_div_delta_M * M) * 1/0.95) + 1

