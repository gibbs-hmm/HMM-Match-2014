%**************************************************************%
% Probabilistic sequence alignment of stratigraphic records    %
%                                                              %
% Please acknowledge the program authors on any publication of %
% scientific results based in part on use of the program and   %
% cite the following articles in which the program was         %
% described.                                                   %
% Lin, L., D. Khider, L. E. Lisiecki and C. E. Lawrence (2014).%
% "Probabilistic sequence alignment of stratigraphic records." %
% Paleoceanography: 2014PA002713.                              %
%                                                              %
% This program is free software; you can redistribute it       %
% and/or modify it under the terms of the GNU General Public   %
% License as published by the Free Software Foundation;        %
% either version 2 of the License, or (at your option)         %
% any later version.                                           %
%                                                              %
% This program is distributed in the hope that it will be      %
% useful, but WITHOUT ANY WARRANTY; without even the implied   %
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR      %
% PURPOSE. See the GNU General Public License for more         %
% details.                                                     %
%                                                              %
% You should have received a copy of the GNU General Public    %
% License along with this program; if not, write to the        %
% Free Software Foundation, Inc., 51 Franklin Street,          %
% Fifth Floor, Boston, MA  02110-1301, USA.                    %
%**************************************************************%

clear

fprintf('1312 original\n')

sequences_info_1312;

%some initial setting up
qqqq=1;
pre_processing_fine;

T = 0;
bb=0;
while T==0  
    tic
    bb=bb+1;
    tic
    forward_point_fine;
    toc
    
    back_sample_point_fine;
    update_parameter_point_script_fmin_fine;
    
    
    fprintf('log likelihood in iteration %d is %f\n',bb, log_fE);
    fprintf('The sigma returned in iteration %d is %f\n', bb, sigma);
    fprintf('The mu returned in iteration %d is %f\n', bb, mu);
    fprintf('The pi_begin returned in iteration %d is %f %f\n', bb, pi_begin(1), pi_begin(2));
    fprintf('The tao_x_begin returned in iteration %d is %f\n', bb, tao_x_begin(1));
    fprintf('The tao_x_end returned in iteration %d is %f\n', bb, tao_x_end(1));
    fprintf('The tao_y_begin returned in iteration %d is %f\n', bb, tao_y_begin(1));
    fprintf('The tao_y_end returned in iteration %d is %f\n', bb, tao_y_end(1));
    fprintf('The delta_r returned in iteration %d is %f %f\n', bb, delta_r(1), delta_r(2));

    determine_stopping;
    toc
end;
clear log_f

%viterbi_point_fine;

uncertainty_anal_point_fine;

output_fine;

clear VV

save(['point_exact_ratio_' num2str(latepl) '.mat'],'-v7.3')
