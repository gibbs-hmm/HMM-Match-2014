fid = fopen('stat.txt','w+');
fprintf(fid,'%s\n',['mean squared error median : ' num2str(scaled_mse)]);
fprintf(fid,'%s\n',['mean CI width : ' num2str(mean(ci_width))]);
fprintf(fid,'%s\n',['median CI width : ' num2str(median(ci_width))]);
fprintf(fid,'%s\n',['Estimated mu is : ' num2str(mu)]);
fprintf(fid,'%s\n',['Estimated sigma is : ' num2str(sigma)]);
fclose(fid);


%plot the raw sequences in target scale
hh = figure;
plot(input_scaled_new(:,1),input_scaled_new(:,2),'r',target(:,1),target(:,2),'b');
legend('input', 'target');
title('sequences before alignment');
xlabel('depth at the site')
ylabel('\delta^{18}O')
saveas(hh,'raw_sequence.eps','epsc')
saveas(hh, 'raw_sequence.fig')


%plot confidence limit for median
confidence_limit_plot(median_ali, 'median', input_track, upper_95, lower_95,latepl,target);

% confidence_limit_plot(cen_ali, 'centroid', input_scaled_new, upper_95, lower_95,latepl);
% 
% confidence_limit_plot(viterbi_ali, 'viterbi', input_scaled_new, upper_95, lower_95,latepl);


%plot aligment and sedimentation rates
[log_prob_med ali_ratio_med log_residual_seq_med ...
    log_rate_change_seq_med rate_seq_med] = ali_seq_anal_plot_gap_new_fine(median_seq,'median',...
    input_scaled_new,target,discrete_target,pi_0,pi_nomatch,...
    dis_bin_input_new,ratio_track, log_grid,mu,latepl,ratio_track_new,transition_prob_track,group_transition);

% [log_prob_cen ali_ratio_cen log_residual_seq_cen ...
%     log_rate_change_seq_cen] = ali_seq_anal_plot_gap_new_fine(cen_seq,'centroid',...
%     input_scaled_new,target,discrete_target,pi_0,pi_nomatch,...
%     dis_bin_input_new,ratio_track, ratio_track_1,log_grid,mu,latepl,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2,group_transition);
% 
% [log_prob_viterbi ali_ratio_viterbi log_residual_seq_viterbi ...
%     log_rate_change_seq_viterbi] = ali_seq_anal_plot_gap_new_fine(viterbi_seq,'viterbi',...
%     input_scaled_new,target,discrete_target,pi_0,pi_nomatch,...
%     dis_bin_input_new,ratio_track, ratio_track_1,log_grid,mu,latepl,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2,group_transition);
