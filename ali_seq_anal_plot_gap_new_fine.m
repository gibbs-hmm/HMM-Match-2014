%given an alignment sequence (contains info of aligned bin index),  plot
%the alignment along with corresponding ratio plot and ratio histogram;
% - filename is description (name) of the alignment
% - also returns a set of information for alignment: 
% - log_prob: log likeligood of the alignment
% - ali_ratio: mapping rate value at each input data point
% - ratio_index_seq: mapping rate index info at each input data point
% - log_residual_seq: log of pdf of residual at each input data point
% - log_rate_change_seq: log of transition probability at each input data
% point

function [log_prob ali_ratio log_residual_seq log_rate_change_seq rate_seq] = ali_seq_anal_plot_gap_new_fine(ali_seq, filename,input_scaled_new,target,discrete_target,pi_0,pi_nomatch,dis_bin_input_new,ratio_track,log_grid, mu,latepl,ratio_track_new,transition_prob_track,group_transition)

[log_prob ali_ratio log_residual_seq log_rate_change_seq] = ...
    calc_ali_prob_fine(ali_seq,input_scaled_new,pi_0,pi_nomatch,dis_bin_input_new,...
    ratio_track,log_grid,ratio_track_new,transition_prob_track,group_transition);
%rmax = length(ratio_r);
ali = zeros(length(input_scaled_new),1);
for i=1:length(input_scaled_new)
    if ali_seq(i)~=0
        ali(i)=discrete_target(ali_seq(i));
    end
end

have_align = find(ali_seq~=0);
align_start=have_align(1);
align_end=have_align(end);

center = floor(mean(target(:,2)));
bottom = floor(min(min(target(:,2)),min(input_scaled_new(:,2))));
if min(min(target(:,2)),min(input_scaled_new(:,2)))>bottom+0.5
    bottom=bottom+0.5;
end
top = ceil(max(max(target(:,2)),max(input_scaled_new(:,2))));
if max(max(target(:,2)),max(input_scaled_new(:,2)))<top-0.5
    top=top-0.5;
end


flip_input = -input_scaled_new(align_start:align_end,2)+2*center;
flip_target = -target(:,2)+2*center;

ytick = -top+2*center:0.5:-bottom+2*center;
ytick_flip = top:-0.5:bottom;



hh=figure;
plot(ali(align_start:align_end),flip_input,'go',target(:,1),flip_target-mu,'b')
legend('alignment','target')
hold on
plot(ali(align_start:align_end),flip_input,'r',target(:,1),flip_target-mu,'b')
legend('alignment','target')
hold off
%title(['scaled alignment plot -- ' filename])
xlabel('Stack Age (kyr)')
ylabel('\delta^{18}O')

set(hh,'Position',[200 200 800 220]);
set(gca, 'YTick', ytick);
set(gca, 'YTickLabel', num2str(ytick_flip'));
axis([0 max(target(:,1)) min(ytick) max(ytick)])

saveas(hh, [filename '_plot_scaled_flip_' num2str(latepl) '.eps'],'epsc')

saveas(hh, [filename '_plot_scaled_flip_' num2str(latepl) '.fig'])

ss = find(ali_ratio~=0);
start_map = ss(1);
tail_map = ss(end);

ali_change=[];
ali_change(1,1)=align_start;
ali_change(1,2)=ali(align_start);
s=1;
old = ali(align_start);
for i=align_start+1:align_end
    if ali(i)~=old
        old = ali(i);
        s=s+1;
        vv = find(ali==old);
        ali_change(s,1)=vv(end);
        ali_change(s,2)=old;
    end
end

rate_seq = zeros(length(ali_change)-1,1);

for i=1:length(rate_seq)
    rate_seq(i)=(input_scaled_new(ali_change(i+1,1),1)-input_scaled_new(ali_change(i,1),1))/...
        (ali_change(i+1,2)-ali_change(i,2));
    pos(i)=ali(ali_change(i+1));
end






hh=figure;
vv=find(ali(2:end)~=0);

% plot(input_scaled(2:end,1),ali_ratio,'b-')
plot(ali(vv),ali_ratio(vv),'b-')
%title(['Relative accumuluatio rate -- ' filename] )
set(hh,'Position',[200 200 800 220]); 
xlabel('Stack Age (kyr)')
ylabel('Accumulation Rates (m/kyr)')
axis([0 max(target(:,1)) 0 4])


saveas(hh, [filename '_ali_ratio_seq_' num2str(latepl) '.eps'],'epsc')

saveas(hh, [filename '_ali_ratio_seq_' num2str(latepl) '.fig'])