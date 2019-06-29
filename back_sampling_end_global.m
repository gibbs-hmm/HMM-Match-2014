
%sampling the end status

n_ali = length(discrete_target)+(length(input_scaled_new)-1)+(length(discrete_target)-1);
log_pp_ali = zeros(n_ali,1);
index_set = zeros(n_ali,2);

vv_ali=0;
vv2=0;
tt = find(gap_add_f_e~=-inf & isnan(gap_add_f_e)==0);
tt2=gap_add_f_e(tt);
tt3=max(tt2);
for i=1:length(tt)
   vv2=vv2+exp(tt2(i)-tt3);
end
vvv=log(vv2)+tt3;

% for i=1:length(discrete_target)
%     if gap_add_f_e(i)~=-inf && isnan(gap_add_f_e(i))==0
%         %vv_ali=vv_ali+1;
%         vv2 =vv2+ exp(gap_add_f_e(i));
%         
%     end
% end

vv_ali=vv_ali+1;
index_set(vv_ali,1:2)=[length(input_scaled_new) length(discrete_target)];
log_pp_ali(vv_ali)=vvv;

sum_gap_add_f_x = zeros(length(input_scaled_new)-1,1);
for i=1:length(input_scaled_new)-1
    tt = max(gap_add_f_x(i,:));
    for j=1:length(discrete_target)
        sum_gap_add_f_x(i) = sum_gap_add_f_x(i) + exp(gap_add_f_x(i,j)-tt);
    end
    sum_gap_add_f_x(i) = log(sum_gap_add_f_x(i))+tt;
    if sum_gap_add_f_x(i)~=-inf && isnan(sum_gap_add_f_x(i))==0
    vv_ali = vv_ali+1;
    log_pp_ali(vv_ali) = sum_gap_add_f_x(i);
    index_set(vv_ali,:)=[i length(discrete_target)];
    end
end

sum_gap_add_f_y = zeros(length(discrete_target)-1,1);
for i=1:length(discrete_target)-1
    tt = max(gap_add_f_y(i,:));
    for j=1:length(discrete_target)-1
        sum_gap_add_f_y(i) = sum_gap_add_f_y(i) + exp(gap_add_f_y(i,j)-tt);
    end
    sum_gap_add_f_y(i) = log(sum_gap_add_f_y(i))+tt;
    if sum_gap_add_f_y(i)~=-inf && isnan(sum_gap_add_f_y(i))==0
    vv_ali = vv_ali+1;
    log_pp_ali(vv_ali) = sum_gap_add_f_y(i);
    index_set(vv_ali,:)=[length(input_scaled_new) i];
    end
end

log_pp_ali = log_pp_ali(1:vv_ali,:);
index_set = index_set(1:vv_ali,:);
end_set_sampling = zeros(1000,3);

for ss=1:1000
    cc = catigorical_random_log(log_pp_ali);
    
    
    ind_i=index_set(cc,1);
    ind_j=index_set(cc,2);
    
    
    if ind_i==length(input_scaled_new) && ind_j==length(discrete_target)
        cc2=catigorical_random_log(gap_add_f_e);
        ind_k=cc2;
    elseif ind_i==length(input_scaled_new)
        cc2 = catigorical_random_log(gap_add_f_y(ind_j,:));
        ind_k=cc2;
    else
        cc2 = catigorical_random_log(gap_add_f_x(ind_i,:));
        ind_k=cc2;
    end
    end_set_sampling(ss,:)=[ind_i ind_j ind_k];
end