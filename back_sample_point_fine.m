


%generate sample alignments -- back sampling
% - sample_seq contains the 1000 samples, each column represents one
% sample,similar for sample_ali
% - sample_seq: records the aligned bin index for alignments
% - sample_ale: records the aligned age info for alignments


sample_ali=zeros(length(input_scaled_new),1000);
sample_seq=zeros(length(input_scaled_new),1000);
%exp_logf=exp(log_f);

back_sampling_end_global;
for nn=1:size(sample_ali,2)
   
    
    curr_i=end_set_sampling(nn,1);
    curr_j=end_set_sampling(nn,2);
    s = end_set_sampling(nn,3);
    sample_seq(curr_i,nn) = curr_j;
    sample_seq(curr_i-1,nn) = s;
    %     tmp = curr_j-s+1;
    %     if dis_bin_input(curr_i)~=0 && tmp~=0
    %     pre_index = ratio_track(dis_bin_input(curr_i),tmp);
    %     elseif tmp~=0
    %         pre_index = ratio_track_0(tmp);
    %
    %     elseif dis_bin_input(curr_i)~=0
    %       pre_index = ratio_track_1(dis_bin_input(curr_i));
    %     else
    %         pre_index = 9;
    %     end
    pre_rr = ratio_track_new(curr_i, curr_j-s+1);
%     if curr_j~=s
%         pre_rr = dis_bin_input_new(curr_i)/(curr_j-s);
%     else
%         pre_rr=dis_bin_input_new(curr_i)/0.5;
%     end
    curr_i=curr_i-1;
    curr_j=s;
    
    
    while curr_i~=1 && curr_j~=1
        u=rand;
        ind_nonzero = find(log_f(curr_i,curr_j,:)~= -Inf);
        
        ssum2=zeros(length(ind_nonzero),1);
        pp=zeros(length(ind_nonzero),1);
        for tt=1:length(ind_nonzero)
            tmp = curr_j-ind_nonzero+1;  % it mapped to
            curr_rr = ratio_track_new(curr_i, ind_nonzero(tt));
%             if curr_j~=ind_nonzero(tt)
%                 curr_rr = dis_bin_input_new(curr_i)/(curr_j-ind_nonzero(tt));
%             else
%                 curr_rr = dis_bin_input_new(curr_i)/0.5;
%             end
            
            %             tmp = curr_j-ind_nonzero(tt)+1;
            %             if dis_bin_input(curr_i)~=0 && tmp~=0
            %             curr_index = ratio_track(dis_bin_input(curr_i),...
            %                 tmp);
            %             elseif tmp~=0
            %                 curr_index = ratio_track_0(tmp);
            %                 %[qq curr_index] = min(abs(ratio_r-dis_bin_input(curr_i))/(curr_j-ind_nonzero(tt)));
            %             elseif dis_bin_input(curr_i)~=0
            %                 curr_index = ratio_track_1(dis_bin_input(curr_i));
            %             else
            %                 curr_index = 9;
            %                % [qq curr_index] = min(abs(ratio_r-dis_bin_input(curr_i)/0.5));
            %             end
            if pre_rr<0.9220
                pre_qq=1;
            elseif pre_rr>1.085
                pre_qq=3;
            else
                pre_qq=2;
            end
            
            if curr_rr<0.9220
                curr_qq=1;
            elseif curr_rr>1.085
                curr_qq=3;
            else
                curr_qq=2;
            end
           % density_rr = density_mixture_gaussian(log(pre_rr),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
          %  bbb1= group_transition(curr_qq,pre_qq)*density_rr;
          density_rr =  transition_prob_track(curr_i,curr_j-s+1);
          bbb1 = group_transition(curr_qq,pre_qq)*density_rr;
            
            if bbb1==0
                pp(tt)=0;
            else
                for t=1:length(ind_nonzero)
%                     if curr_j~=ind_nonzero(t)
%                         curr_rr2=dis_bin_input_new(curr_i)/(curr_j-ind_nonzero(t));
%                     else
%                         curr_rr2=dis_bin_input_new(curr_i)/0.5;
%                     end
                    curr_rr2 = ratio_track_new(curr_i,ind_nonzero(t));
                   
                    %                     tmp = curr_j-ind_nonzero(t)+1;
                    %                     if dis_bin_input(curr_i)~=0 && tmp~=0
                    %                         curr_index2 = ratio_track(dis_bin_input(curr_i),tmp);
                    %                     elseif tmp~=0
                    %                         curr_index2 = ratio_track_0(tmp);
                    %                         %[qq curr_index2] = min(abs(ratio_r-dis_bin_input(curr_i))/(curr_j-ind_nonzero(t)));
                    %                     elseif dis_bin_input(curr_i)~=0
                    %                         curr_index2 = ratio_track_1(dis_bin_input(curr_i));
                    %                     else
                    %                         curr_index2 = 9;
                    %                         %[qq curr_index2] = min(abs(ratio_r-dis_bin_input(curr_i)/0.5));
                    %                     end
                    if curr_rr2<0.9220
                        curr_qq2=1;
                    elseif curr_rr2>1.085
                        curr_qq2=3;
                    else
                        curr_qq2=2;
                    end
                    
                    bbb2= group_transition(curr_qq2,pre_qq)*density_rr;
                    
                    
                    if bbb2~=0
                        ssum2(tt)=ssum2(tt)+exp(log_f(curr_i,curr_j,...
                            ind_nonzero(t))-log_f(curr_i,curr_j,ind_nonzero(tt)))...
                            *bbb2/bbb1;
                    end;
                end;
                pp(tt)=1/ssum2(tt);
            end;
        end;
        
        candidate=find(pp~=0);
        ssss=0;
        for s=1:size(candidate,1)
            ssss=ssss+pp(candidate(s));
            if u<=ssss
                %  sample_seq(curr_i-1,nn)=
                sample_seq(curr_i-1,nn)=curr_j-ind_nonzero(candidate(s))+1;
                break
                
            end
        end
        s=curr_j-ind_nonzero(candidate(s))+1;
        
        pre_rr = ratio_track_new(curr_i,curr_j-s+1);
%         if curr_j~=s
%             pre_rr = dis_bin_input_new(curr_i)/(curr_j-s);
%         else
%             pre_rr = dis_bin_input_new(curr_i)/0.5;
%         end
%         
        %         tmp = curr_j-s+1;
        %         if dis_bin_input_new(curr_i)~=0 && tmp~=0
        %             pre_index = ratio_track(dis_bin_input_new(curr_i),tmp);
        %         elseif tmp~=0
        %             pre_index = ratio_track_0(tmp);
        %             %[qq pre_index] = min(abs(ratio_r-dis_bin_input(curr_i)/(curr_j-s)));
        %         elseif dis_bin_input(curr_i)~=0
        %             pre_index = ratio_track_1(dis_bin_input(curr_i));
        %         else
        %             pre_index = 9;
        %             %[qq pre_index] = min(abs(ratio_r-dis_bin_input(curr_i)/0.5));
        %         end
        curr_i=curr_i-1;
        curr_j=s;
        
        
        
        
    end
    
    
end

for nn=1:size(sample_ali,2)
    for i=1:length(input_scaled_new)
        if sample_seq(i,nn)>0
            sample_ali(i,nn)=discrete_target(sample_seq(i,nn));
        end
    end
end