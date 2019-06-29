function residual_seq = calc_ali_residual_fine(ali_seq,input_scaled_new,dis_bin_input_new,delta_grid,ratio_track_new)

ali_ratio = zeros(length(input_scaled_new)-1,1);
residual_seq = zeros(length(input_scaled_new),1);


have_align = find(ali_seq~=0);
align_start = have_align(1);
align_end = have_align(end);
% if align_start>1
%     for i=1:align_start-1
%     residual_seq(i)=input(i,2)-mean(input(:,2));
%     end
% end
% 
% if align_end<length(input)
%     
%     for i=align_end+1:length(input)
%         residual_seq(i)=input(i,2)-mean(input(:,2));
%     end
% end
    
    
residual_seq(align_start) = delta_grid(align_start,ali_seq(align_start));
i = align_start+1;
% tmp = ali_seq(i)-ali_seq(i-1);
% if tmp~=0 && dis_bin_input_new(i)~=0
%     pre_index = ratio_track(dis_bin_input_new(i),tmp);
% %elseif tmp~=0
%  %   pre_index = ratio_track_0(tmp);
% elseif dis_bin_input_new(i)~=0
%     pre_index = ratio_track_1(dis_bin_input_new(i));
% else
%     pre_index = 9;
% end

% if ali_seq(i)~=ali_seq(i-1)
%     pre_rr = dis_bin_input_new(i)/(ali_seq(i)-ali_seq(i-1));
% else
%     pre_rr = dis_bin_input_new(i)/0.5;
% end

pre_rr = ratio_track_new(i,ali_seq(i)-ali_seq(i-1)+1);
    
%pre_index = ratio_track(dis_bin_input_new(i),ali_seq(i)-ali_seq(i-1)+1);
%ratio_index_seq(i-1) = pre_index;
ali_ratio(i-1)=pre_rr;

residual_seq(i)=delta_grid(i,ali_seq(i));

for i=align_start+2:align_end
%     tmp = ali_seq(i)-ali_seq(i-1);
%     if tmp~=0 && dis_bin_input_new(i)~=0
%         curr_index = ratio_track(dis_bin_input_new(i),tmp);
%     elseif tmp~=0
%         curr_index = ratio_track_0(tmp);
%     elseif dis_bin_input_new(i)~=0
%         curr_index = ratio_track_1(dis_bin_input_new(i));
%         
%     else
%         curr_index = 9;
%     end
%     
%     if ali_seq(i)-ali_seq(i-1)~=0
%         curr_rr = dis_bin_input_new(i)/(ali_seq(i)-ali_seq(i-1));
%     else
%         curr_rr = dis_bin_input_new(i)/0.5;
%     end

curr_rr = ratio_track_new(i,ali_seq(i)-ali_seq(i-1)+1);
    %curr_index = ratio_track(dis_bin_input_new(i),ali_seq(i)-ali_seq(i-1)+1);
    
    pre_rr = curr_rr;
    ali_ratio(i-1) = pre_rr;
    residual_seq(i)=delta_grid(i,ali_seq(i));
    
end

% ali_ratio = zeros(length(ratio_index_seq),1);
% for i=1:length(ratio_index_seq)
%     if ratio_index_seq(i)~=0
%         ali_ratio(i)=ratio_r(ratio_index_seq(i));
%     end
% end