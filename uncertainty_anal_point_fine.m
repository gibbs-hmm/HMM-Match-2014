%obtain median estimator

calc_median_ali_dp_fine_half;

%calculate upper and lower limit
[upper_95 lower_95] = upper_lower_limit(input_scaled_new, sample_ali);


residual_seq_median = calc_ali_residual_fine(median_seq,input_scaled_new,dis_bin_input_new,delta_grid,ratio_track_new);
ci_width = upper_95 - lower_95;

vv = find(residual_seq_median==0);
tt = residual_seq_median;
tt(vv)=input_scaled_new(vv,2)-mean(input_scaled_new(:,2));

tt=residual_seq_median(find(residual_seq_median~=0));
scaled_mse = mean((tt-mu).^2);

