%% Roi-based MVPA for single subject (run_split_half_correlations_single_sub)
%
% Load t-stat data from one subject, apply 'vt' mask, compute difference
% of (fisher-transformed) between on- and off diagonal split-half
% correlation values.
%
% #   For CoSMoMVPA's copyright information and license terms,   #
% #   see the COPYING file distributed with CoSMoMVPA.           #

clear all
addpath S:\nad12\CoSMoMVPA-master 
%% Set analysis parameters
subject_id='18y404'; %
roi_label='LeftMTL'; % name of ROI mask used for running correlations

%config=cosmo_config(); test on and off
study_path='S:\nad12\FAME8\Analysis_ret\FAME_categorymodel_ret_hrf';

%% Computations
data_path=fullfile('S:\nad12\FAME8\Analysis_ret\FAME_categorymodel_ret_hrf', subject_id);
%sepfiles to combine

% file locations for both halves
% for categort model this will be targets and related lures and also
% unrelatedd lures as a control condition of sorts....
%so I am guesisng we will need to adjust the 1 and 2 numbers below, right
%Indira? else we just want to look at 2 trials and not 3
%  ** AvgTargetBetas.nii; AvgRelLureBetas.nii; AvgUnrelLureBetas.nii  **
Enc_Rem=fullfile(data_path, 'Enc_AvgRemBetas.nii');
Enc_Know=fullfile(data_path, 'Enc_AvgKnowBetas.nii');
Ret_Rem=fullfile(data_path, 'Ret_AvgRemBetas.nii');
Ret_Know=fullfile(data_path, 'Ret_AvgKnowBetas.nii');
mask_fn=fullfile(data_path,[roi_label '_Mask.nii']); %second half of mask name

% load two halves as CoSMoMVPA dataset structs.
%Chunks = Runs  Targets = trial type conditions
Enc_Rem_ds=cosmo_fmri_dataset(Enc_Rem,'mask',mask_fn,... %encoding... run 1 (chunk) Target(rem) will be identified as 1
                                     'targets',(1)',... %ex: Rem
                                     'chunks',(1)); %ex: encoding (avg. of all enc runs)
                                 
Enc_Know_ds=cosmo_fmri_dataset(Enc_Know,'mask',mask_fn,...
                                     'targets',(2)',... %ex. Know
                                     'chunks',(1)); %encoding (avg. of all enc runs)
                                 
Ret_Rem_ds=cosmo_fmri_dataset(Ret_Rem,'mask',mask_fn,...
                                     'targets',(1)',... %Rem
                                     'chunks',(2)); %ret
                                 
Ret_Know_ds=cosmo_fmri_dataset(Ret_Know,'mask',mask_fn,...
                                     'targets',(2)',... %Know
                                     'chunks',(2)); %ret

%Combine files at encoding and retrieval to create two files (i.e.,
%stacking)
%make sure all ds_* changed from here on
ds_enc=cosmo_stack({Enc_Rem_ds, Enc_Know_ds}); 
ds_ret=cosmo_stack({Ret_Rem_ds, Ret_Know_ds});

labels_enc={'Rem_enc'; 'Know_enc'}; %outcome of condition labeling
labels_ret={'Rem_ret'; 'Know_ret'};
ds_enc.sa.labels = labels_enc;
ds_ret.sa.labels = labels_ret;

%cosmo fxn to make sure data in right format
cosmo_check_dataset(ds_enc);
cosmo_check_dataset(ds_ret);

% Some sanity checks to ensure that the data has matching features (voxels)
% and matching targets (conditions)
assert(isequal(ds_enc.fa,ds_ret.fa));
assert(isequal(ds_enc.sa.targets,ds_ret.sa.targets));

%change if you change ds_* above
nClasses = numel(ds_enc.sa.labels);

% get the sample data - samples are the correlations being ran on each
% voxel (a bunch of numbers)
ds_enc_samples=ds_enc.samples;
ds_ret_samples=ds_ret.samples;

% compute all correlation values between the two halves, resulting
% in a 6x6 matrix. Store this matrix in a variable 'rho'.
% Hint: use cosmo_corr (or builtin corr, if the matlab stats toolbox
%       is available) after transposing the samples in the two halves.
% >@@>
rho=cosmo_corr(ds_enc_samples',ds_ret_samples');
% <@@<

% Correlations are limited between -1 and +1, thus they cannot be normally
% distributed. To make these correlations more 'normal', apply a Fisher
% transformation and store this in a variable 'z'
% (hint: use atanh).
% >@@>
z=atanh(rho);
% <@@<

%% Save Workspace, Clear Workspace, Redefine subjects cell array.

save('S:\nad12\FAME8\Analysis_ret\FAME_categorymodel_ret_hrf\ERS_Targ_Lure_Workspace.mat')

%% visualize the normalized correlation matrix
figure
% >@@>
imagesc(z);
colorbar()
set(gca, 'xtick', 1:numel(ds_enc.sa.labels), ...
                'xticklabel', ds_enc.sa.labels)
set(gca, 'ytick', 1:numel(ds_ret.sa.labels), ...
                'yticklabel', ds_ret.sa.labels)
title(subject_id)
% <@@<

% Set up a contrast matrix to test whether the element in the diagonal
% (i.e. a within category correlation) is higher than the average of all
% other elements in the same row (i.e. the average between-category
% correlations). For testing the split half correlation of n classes one
% has an n x n matrix (here, n=6).
%
% To compute the difference between the average of the on-diagonal and the
% average of the off-diagonal elements, consider that there are
% n on-diagonal elements and n*(n-1) off-diagonal elements.
% Therefore, set
% - the on-diagonal elements to 1/n           [positive]
% - the off-diagonal elements to -1/(n*(n-1)) [negative]
% This results in a contrast matrix with weights for each element in
% the correlation matrix, with positive and equal values on the diagonal,
% negative and equal values off the diagonal, and a mean value of zero.
%
% Under the null hypothesis one would expect no difference between the
% average on the on- and off-diagonal, hence correlations weighted by the
% contrast matrix has an expected mean of zero. A postive value for
% the weighted correlations would indicate more similar patterns for
% patterns in the same condition (across the two halves) than in different
% conditions.

% Set the contrast matrix as described above and assign it to a variable
% named 'contrast_matrix'
% >@@>
contrast_matrix=(eye(nClasses)-1/nClasses)/(nClasses-1);

% alternative solution
contrast_matrix_alt=zeros(nClasses,nClasses);
for k=1:nClasses
    for j=1:nClasses
        if k==j
            value=1/nClasses;
        else
            value=-1/(nClasses*(nClasses-1));
        end
        contrast_matrix_alt(k,j)=value;
    end
end

% <@@<

% sanity check: ensure the matrix has a sum of zero
if abs(sum(contrast_matrix(:)))>1e-14
    error('illegal contrast matrix: it must have a sum of zero');
end

%visualize the contrast matrix
% >@@>
figure
imagesc(contrast_matrix)
colorbar
title('Contrast Matrix')
% <@@<

% Weigh the values in the matrix 'z' by those in the contrast_matrix
% and then average them (hint: use the '.*' operator for element-wise
% multiplication).
% Store the result in a variable 'weighted_z'.
% >@@>
weighted_z=z .* contrast_matrix;
% <@@<

% Compute the sum of all values in 'weighted_z', and store the result in
% 'sum_weighted_z'.
% >@@>
sum_weighted_z=sum(weighted_z(:)); %Expected value under H0 is 0
% <@@<

%visualize weighted normalized correlation matrix
figure
imagesc(weighted_z)
colorbar
% For the advanced exercise
title(sprintf('Weighted Contrast Matrix m = %5.3f', sum_weighted_z))

%Create code to output file... 
%Look at weighted Z