diff = abs(mu_old-mu);
diff = max(diff,abs(sigma-sigma_old));

diff2 = max(abs(tao_x_begin-tao_x_begin_old));
diff = max(diff,diff2);

diff = max(diff, max(abs(tao_x_begin-tao_x_begin_old)));
diff = max(diff, max(abs(tao_x_end-tao_x_end_old)));
diff = max(diff, max(abs(tao_y_begin-tao_y_begin_old)));
diff = max(diff, max(abs(tao_y_end-tao_y_end_old)));

diff2 = max(abs(delta_r-delta_r_old));
diff = max(diff2,diff);

diff2 = max(abs(pi_begin_old - pi_begin));
diff = max(diff2,diff);

if diff<0.001 || bb>10
    T=1;
else
    T=0;
end

mu_old = mu;
sigma_old = sigma;
tao_x_begin_old = tao_x_begin;
tao_x_end_old = tao_x_end;
tao_y_begin_old = tao_y_begin;
tao_y_end_old = tao_y_end;
delta_r_old = delta_r;
pi_begin_old = pi_begin;