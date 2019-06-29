%updated the way to calculate the rates


%preprocessing the input and target sequences
%select the sequence from the whole record
seq1 = ori1(min(find(ori1(:,1)>=begin1)):max(find(ori1(:,1)<=end1)),:);
seq2 = ori2(min(find(ori2(:,1)>=begin2)):max(find(ori2(:,1)<=end2)),:);

%filter out duplicate data
[input(:,1),label]=unique(seq1(:,1));
input(:,2)=seq1(label,2);
target=seq2;

%scale the starting point to zero
target_start=target(1,1);
target(:,1)=target(:,1)-target_start;

%scale the input such that it has the same scale as target
input_scaled = input;
input_scaled(:,1)=target(1,1)+(input(:,1)-input(1,1))*(target(end,1)-target(1,1))/...
    (input(end,1)-input(1,1));

tt= load('sedrate_dist_evenbins.txt');
marginal_dis = load('prob_original');
group_count = zeros(3,3);
group_count(1,1)=114;
group_count(1,2)=21;
group_count(1,3)=27;
group_count(2,1)=14;
group_count(2,2)=23;
group_count(2,3)=25;
group_count(3,1)=29;
group_count(3,2)=18;
group_count(3,3)=184;



   
   
%record the ratio values for each ratio
ratio_r = tt(:,2);
min_r = 0.25;
max_r = 4;

%ratio=[1,4;1,3;2,5;1,2;3,5;2,3;3,4;4,5;1,1;5,4;4,3;3,2;5,3;2,1;5,2;3,1;4,1];
rmax = length(ratio_r);

mix_std1 = sqrt(0.0216);
mix_std2 = sqrt(0.0929);
mix_p1 = 0.642432;
mix_p2 = 1-mix_p1;
mix_mu1 = 0.0198;
mix_mu2 = -0.0297;

%0.642432*normpdf(q,0.0198, sqrt(0.0216))+0.357568*normpdf(q,-0.0297,sqrt(0.0929))
%  density_mixture_gaussian(x,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2)

 tt = linspace(log10(0.25),log10(0.9220),50);
   s=0;
   hh=tt(2)-tt(1);
   for j=1:length(tt)-1
       q=(tt(j)+tt(j+1))/2;
       %kk1=0.642432*normpdf(q,0.0198, sqrt(0.0216))+0.357568*normpdf(q,-0.0297,sqrt(0.0929));
       kk1=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
       s=s+abs(kk1);
   end
   left_marginal=s*hh;
   
    tt = linspace(log10(1.0850),log10(4),100);
   s=0;
   hh=tt(2)-tt(1);
   for j=1:length(tt)-1
       q=(tt(j)+tt(j+1))/2;
       %kk1=0.642432*normpdf(q,0.0198, sqrt(0.0216))+0.357568*normpdf(q,-0.0297,sqrt(0.0929));
       kk1=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
       s=s+abs(kk1);
   end
   right_marginal=s*hh;
   
    tt = linspace(log10(0.9220),log10(1.0850),100);
   s=0;
   hh=tt(2)-tt(1);
   for j=1:length(tt)-1
       q=(tt(j)+tt(j+1))/2;
       %kk1=0.642432*normpdf(q,0.0198, sqrt(0.0216))+0.357568*normpdf(q,-0.0297,sqrt(0.0929));
       kk1=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
       s=s+abs(kk1);
   end
  mid_marginal=s*hh;
   
   
   log_left_marginal = log(left_marginal);
   log_right_marginal = log(right_marginal);
   log_mid_marginal = log(mid_marginal);
   

sigma = 0.25;
mu = 0;

pi_nomatch = 0;



n1=size(input_scaled,1);
n2=size(target,1);
beginpoint=target(1);
endpoint=target(n2);

%hh is the minimum distance between input points
hh=inf;
dist = zeros(length(input_scaled)-1,1);
for i=2:length(input_scaled)
    dist(i-1)=input_scaled(i,1)-input_scaled(i-1,1);
    if input_scaled(i,1)-input_scaled(i-1,1)<hh
        hh=input_scaled(i,1)-input_scaled(i-1,1);
    end
end

[kk1 kk2] = sort(dist);


%select bin width sunch that each bin contains at least one data point
h3 = target(2,1)-target(1,1);
for i=3:length(target)
    if target(i,1)-target(i-1,1)>h3
        h3 = target(i,1)-target(i-1,1);
    end
end



h = kk1(floor(length(kk1)*0.1))/4;
num_tar = target(end,1)/h;
%h = target(end,1)/ceil(num_tar);

discrete_target = linspace(target(1,1),target(end,1),(ceil(num_tar)+1)/qqqq)';
h = discrete_target(2)-discrete_target(1);

gap_con = ceil(0.02*length(input));

gap_con_tar = ceil(10*h3/h);


% discrete_target=[];
% discrete_target(1,1)=target(1,1);
% ss=1;
% while discrete_target(ss,1)<target(end,1)
%    ss=ss+1;
%    discrete_target(ss,1)=discrete_target(ss-1,1)+h;
% end


bin_mean = zeros(length(discrete_target),1);
bin_mean(1)=target(1,2);
for i=2:length(target)
   tmp = find(discrete_target>target(i-1,1) & discrete_target<=target(i,1));
   for j=1:length(tmp)
      bin_mean(tmp(j))=target(i-1,2)+(target(i,2)-target(i-1,2))*...
          (discrete_target(tmp(j))-target(i-1,1))/(target(i,1)-target(i-1,1)); 
   end
end


% %discretize the target into corresponding bins
% discrete_target = (target(:,1):h:target(end,1))';
% interval_target_size=zeros(length(discrete_target),1);  %number of points in each bin
% interval_target = cell(length(discrete_target),1);  %record the index for points in each bin
% a=0;
% for i=1:length(interval_target_size)
%     temp=target(a+1:size(target,1),1);
%     ttt=find(temp<beginpoint+i*h)+a;
%     if length(ttt)~=0
%         a=max(ttt);
%     end;
%     interval_target{i}=ttt;
%     interval_target_size(i)=length(ttt);
%     clear ttt
% end;

%assign the index of bin for input when putting input in discreted target bins
input_bin_index = zeros(length(input),1);
ss =1;
for i = 1: length(discrete_target)
    while input_scaled(ss,1)<=discrete_target(i)
        input_bin_index(ss)=i;
        ss = ss+1;
        if ss>length(input)
            break
        end
    end
end
input_bin_index(ss:end)=length(discrete_target);


%record the number of bin distance betwen input point with previous point
%set the distance for the first element as 0
dis_bin_input = zeros(length(input),1);
for i=2:length(dis_bin_input)
    dis_bin_input(i) = input_bin_index(i) - input_bin_index(i-1);
end



input_scaled_new = input_scaled;
input_bin_index_new = input_bin_index;
dis_bin_input_new = dis_bin_input;
tt = find(dis_bin_input==0);
tt = setdiff(tt,[1]);
input_track = input(:,1);
input_track_2 = input(:,2);

for i=1:length(tt)
    input_scaled_new(tt(i)-1,1:2) = (input_scaled(tt(i),1:2)+input_scaled(tt(i)-1,1:2))/2;
    input_track(tt(i)-1)=(input(tt(i),1)+input(tt(i)-1,1))/2;
    input_track_2(tt(i)-1)=(input(tt(i),2)+input(tt(i)-1,2))/2;
end

tmp1 = input_track(:,1);

tmp1(tt)=[];
input_track = tmp1;


tmp1 = input_scaled_new(:,1);
tmp2=input_scaled_new(:,2);
tmp1(tt)=[];
tmp2(tt)=[];


tmp2=input_track_2(:,1);
tmp2(tt)=[];
input_track_2=tmp2;


input_bin_index_new(tt)=[];
dis_bin_input_new(tt)=[];
% for i=1:length(tt)
%    tmp1(tt(i))=[];
%    tmp2(tt(i))=[];
%    input_bin_index_new(tt(i))=[];
%    dis_bin_input_new(tt(i))=[];
% end

dis_bin_input_new = dis_bin_input_new+1;
input_scaled_new = zeros(length(tmp1),2);
input_scaled_new(:,1)=tmp1;
input_scaled_new(:,2)=tmp2;


% 
% %calculate the mean value for each bin
% bin_mean = zeros(length(discrete_target),1);
% for i=1:length(discrete_target)
%     ll1 = interval_target{i}(1);
%     ll2 = interval_target{i}(end);
%     bin_mean(i) = sum(target(ll1:ll2,2))/(ll2-ll1+1);
% end
% 





%look up table -- record the log(pdf) and residue if data point i in input
%is aligned to j^th bin in target
log_grid = zeros(length(input_scaled_new),length(bin_mean));
delta_grid = zeros(length(input_scaled_new),length(bin_mean));
for i=1:length(input_scaled_new)
    for j=1:length(bin_mean)
        %log_grid(i,j) = log(normpdf(input_scaled_new(i,2)-bin_mean(j),mu,sigma));
        delta_grid(i,j) = input_scaled_new(i,2)-bin_mean(j);
    end
end

log_grid = log(normpdf(delta_grid, mu,sigma));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting the starting parameters



%setting uniform distribution as starting for EM algorithm
pi_0 = ones(rmax,1)/rmax;



group_dis = zeros(3,3);
for i=1:3
    group_dis(i,:) = group_count(i,:)/sum(group_count(i,:));
end


A_prior = zeros(rmax,rmax);
g11 = group_dis(1,1)*marginal_dis(1:8)/sum(marginal_dis(1:8));
g13 = group_dis(1,3)*marginal_dis(10:17)/sum(marginal_dis(10:17));
for i=1:8
    A_prior(i,1:8) = g11;
    A_prior(i,9) = group_dis(1,2);
    A_prior(i,10:17) = g13;
end
A_prior(9,1:8) = group_dis(2,1)*marginal_dis(1:8)/sum(marginal_dis(1:8));
A_prior(9,9) = group_dis(2,2);
A_prior(9,10:17) = group_dis(2,3)*marginal_dis(10:17)/sum(marginal_dis(10:17));

g31 = group_dis(3,1)*marginal_dis(1:8)/sum(marginal_dis(1:8));
g33 = group_dis(3,3)*marginal_dis(10:17)/sum(marginal_dis(10:17));
for i=10:17
    A_prior(i,1:8) = g31;
    A_prior(i,9) = group_dis(3,2);
    A_prior(i,10:17) = g33;
end
   
A = A_prior;
log_A = log(A);

%transitions from begin node to x, rate, y
pi_begin = [0.05 0.9 0.05];

log_pi_begin = log(pi_begin);

%transitions from x^begin , y^begin, y^end, x^end
tao_x_begin = [0.1 0.9];
tao_y_begin = [0.1 0.9];

tao_x_end = [0.1 0.9];
tao_y_end = [0.1 0.9];
log_tao_x_begin = log(tao_x_begin);
log_tao_y_begin = log(tao_y_begin);
log_tao_x_end = log(tao_x_end);
log_tao_y_end = log(tao_y_end);
%transitions from the rate state
delta_r = [0.9 0.05 0.05];
log_delta_r = log(delta_r);


mean_tar=mean(target(:,2));
log_input_gap_raw = zeros(length(input_scaled_new),1);
for i=1:length(input_scaled_new)
    log_input_gap_raw(i) = log(normpdf(input_scaled_new(i,2)-mean_tar,mu,sigma));
end


%record starting parameters as old
mu_old = mu;
sigma_old = sigma;
pi_begin_old = pi_begin;
tao_x_begin_old = tao_x_begin;
tao_x_end_old = tao_x_end;
tao_y_begin_old = tao_y_begin;
tao_y_end_old = tao_y_end;
delta_r_old = delta_r;





log_transition_prob_track = zeros(length(input_scaled_new),max(dis_bin_input_new)*4);

transition_prob_track = zeros(length(input_scaled_new),max(dis_bin_input_new)*4);
for i=2:size(transition_prob_track,1)
    for j=1:size(transition_prob_track,2)
       
            rr = dis_bin_input_new(i)/(j);
            log_rr = log10(rr);
            transition_prob_track(i,j)=density_mixture_gaussian(log_rr,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
            log_transition_prob_track(i,j)=log(transition_prob_track(i,j));
       
    end
end

log_group_transition = log(group_dis);
log_group_transition(:,1)=log_group_transition(:,1)-log_left_marginal;
log_group_transition(:,2)=log_group_transition(:,2)-log_mid_marginal;
log_group_transition(:,3)=log_group_transition(:,3)-log_right_marginal;


group_transition = exp(log_group_transition);

ratio_track = zeros(max(dis_bin_input_new),4*max(dis_bin_input_new));
for i=1:max(dis_bin_input_new)
    for j=1:4*i
        [aa ratio_track(i,j)] = min(abs(ratio_r-(i-1)/max(j-1,0.5)));
    end
end

ratio_track_1=zeros(max(dis_bin_input_new),1);
for i=1:max(dis_bin_input_new)
    %if i/0.5>=min(ratio_r) && i/0.5<=max(ratio_r)
    [aa ratio_track_1(i)]=min(abs(ratio_r-i/0.5));
    %end
end

ratio_track_new = zeros(length(input_scaled_new),4*max(dis_bin_input_new));
for i=2:length(input_scaled_new)
    for j=1:4*max(dis_bin_input_new)
        if j~=1
        ratio_track_new(i,j) = (dis_bin_input_new(i)-1)/(j-1);
        else
           ratio_track_new(i,j)=(dis_bin_input_new(i)-1)/0.5; 
        end
    end
end
