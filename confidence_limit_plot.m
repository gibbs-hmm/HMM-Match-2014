%plot confidence limit for given alignment
% - filename is the description (or name) of the alignment
% - ali: alignment with aligned age info
function confidence_limit_plot(ali, filename, input, upper_95, lower_95,latepl,target) 


aa = find(ali~=0);
align_start=aa(1);
align_end=aa(end);

ali_new=ali;
for i=align_end+1:length(ali)
    ali_new(i)=ali(align_end);
end

hh=figure;
plot(ali(align_start:align_end),input(align_start:align_end,1),'r',upper_95,input(:,1),'b',lower_95,input(:,1),'g')
legend(filename,'upper limit', 'lower limit')
%title('point-version -- 95% confidence limit of aligned age at input')
xlabel('Estimated Age (kyr)')
ylabel('Depth in input (m)')

tt1=max(ali(align_start:align_end));
tt2=max(upper_95);
tt3=max(lower_95);
axis([0 max([tt1 tt2 tt3]) 0 max(input(:,1))])

axis([0 max(target(:,1)) 0 max(input(:,1))])
%set(hh,'Position',[200 200 800 220]);
% saveas(hh,['align_age_95_band_' filename '_' num2str(latepl) '.eps'],'epsc')
% saveas(hh, ['align_age_95_band_' filename '_' num2str(latepl) '.fig'])
% 

hh=figure;
plot(ali_new,ali-ali,'r',ali_new,upper_95-ali_new,'b',ali_new,...
    lower_95-ali_new,'g')
legend(filename,'upper limit', 'lower limit')
%title('point-version -- 95% confidence limit of aligned age (normalization) at input')
xlabel('Estimated Age (kyr)')
ylabel('Confidence Band in Age (kyr)')

%axis([0 max([tt1 tt2 tt3]) -5 10]);
%axis([0 max(target(:,1)) -15 10]);
set(hh,'Position',[200 200 800 220]);
% saveas(hh,['align_age_95_band_norm_' filename '_' num2str(latepl) '.eps'],'epsc')
% saveas(hh, ['align_age_95_band_norm_' filename '_' num2str(latepl) '.fig'])