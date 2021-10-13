%%
%Author: Marisa C. Ross, PhD
%Date: March 2021


clear

%%
%this script calls other relevant functions necessary for doing large-scale
%network analysis using asymptotical surprise and the PACO algorithm:
%%%%%%%%%%%%%%%
%1) identify ROIs that are common to all participants; ie., remove ROIs that
%are outside the brain / not covered in some participants
%%%%%%%%%%%%%%%
%2) stack subject correlation matrices together and calculate sample-level
%correction matrix, using ComBat to harmonize across sites (i.e. big_r)
%%%%%%%%%%%%%%%
%3) Do percolation analysis to properly threshold correlation matrix and remove
%negative connections
%Apply PACO algorithm to calculate sample-level community structure and
%generate sample-level spatial map
%%%%%%%%%%%%%%%
%4) do within subject network analyses using the group-level mask, calculating large-scale graph
%theory indices and storing connectivity matrices and subject specific
%spatial maps of the resulting statistics
%%%%%%%%%%%%%%%

%%

sub_ids = load('path/to/subject/file');

%specify location of atlas
atlas_location='/path/to/atlas/file';

%create cell array of each participant's path location and filenames to load
%into function later
for sn=1:numel(sub_ids)
    
    if sub_ids(sn)<=41 %Tours
        if sub_ids(sn)<10
            sub_num=['000' num2str(sub_ids(sn))];
        elseif sub_ids(sn)>=10&&sub_ids(sn)<100
            sub_num=['00' num2str(sub_ids(sn))];
        end
        resting_filename(sn,1)={['/path/to/site/directory/' sub_num '/preproc.scaled.resid+orig']};
    
    elseif sub_ids(sn) >= 42 && sub_ids(sn) <= 126 %Emory
        if sub_ids(sn)>=10&&sub_ids(sn)<100
            sub_num=['00' num2str(sub_ids(sn))];
        elseif sub_ids(sn) >=100 
            sub_num=['0' num2str(sub_ids(sn))];
        end
        
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']};  
    
    elseif sub_ids(sn)  >= 127 && sub_ids(sn) <=206 %McLean 
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory/' sub_num '/preproc.scaled.resid+orig']};
   
    elseif sub_ids(sn) >= 207 && sub_ids(sn) <=256 %Vanderbilt
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory/' sub_num '/preproc.scaled.resid+orig']};
    
    elseif sub_ids(sn) >= 257 && sub_ids(sn) <=362 %UWCisler
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory/' sub_num '/preproc.scaled.resid+orig']};

      elseif sub_ids(sn) >= 363 && sub_ids(sn) <=395 %Gronigen
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory/' sub_num '/preproc.scaled.resid+orig']};          

    elseif sub_ids(sn) >= 445 && sub_ids(sn) <=507 %Michigan
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']};       
                
    elseif sub_ids(sn) >= 508 && sub_ids(sn) <=551 %Munster
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']};    
    
    elseif sub_ids(sn) >= 552 && sub_ids(sn) <=633 %Duke
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']};       
    
    elseif sub_ids(sn) >= 634 && sub_ids(sn) <=673 %UMN
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']};       

    elseif sub_ids(sn) >= 674 && sub_ids(sn) <=731 %UWGrupe
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']};       
    
    elseif sub_ids(sn) >= 732 && sub_ids(sn) <=790 %AMC
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']};       
    
    elseif sub_ids(sn) >= 791 && sub_ids(sn) <=868 %Milwaukee
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    
    elseif sub_ids(sn) >= 869 && sub_ids(sn) <=927 %Utrecht
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    
    elseif sub_ids(sn) >= 928 && sub_ids(sn) <=999 %Minneapolis VA
        sub_num=['0' num2str(sub_ids(sn))];
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
       
    elseif sub_ids(sn) >= 1000 && sub_ids(sn) <=1070 %Stanford Brains
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
        
    elseif sub_ids(sn) >= 1071 && sub_ids(sn) <=1211 %Western Ontario
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    
    elseif sub_ids(sn) >= 1212 && sub_ids(sn) <=1271  %Stanford CausCon
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    elseif sub_ids(sn) >= 1271 && sub_ids(sn) <=1334  %Masaryk
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    elseif sub_ids(sn) >= 1335 && sub_ids(sn) <=1408  %Columbia
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    elseif sub_ids(sn) >= 1409 && sub_ids(sn) <=1459  %Ghent
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    elseif sub_ids(sn) >= 1460 && sub_ids(sn) <= 1499  %Toledo MVA
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    elseif sub_ids(sn) >= 1500 && sub_ids(sn) <=1526  %Toledo ONG
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
    elseif sub_ids(sn) >= 1537 %Waco
        sub_num=num2str(sub_ids(sn));
        resting_filename(sn,1)={['/path/to/site/directory' sub_num '/preproc.scaled.resid+orig']}; 
       
    end 
    
end

%% call first function to identify common ROIs across resting-state scans
tic
result_step_one_surprise_250 =find_common_ROIs_250_function(sub_ids,atlas_location,resting_filename);
toc
%save output of this function for later use, as all subsequent scripts will
%use it
%this takes about 81 minutes with 836 subs
save result_step_one_surprise result_step_one_surprise_250

good_ROIs = result_step_one_surprise_250.good_ROIs;
mask_values = result_step_one_surprise_250.mask_values;
%step a here: stack subject matrices and caculate big r
%choose covariates of interest
clinical_data=load(['/path/to/clinical_data']);
keep_clinical=find(ismember(clinical_data(:,1),sub_ids)==1);
clinical_data=clinical_data(keep_clinical,:);

age = clinical_data(:,3); 
ptsd_severity=clinical_data(:,7); 
sex = clinical_data(:,4);
ptsd_dx = clinical_data(:,5);
group = clinical_data(:,6);
index_trauma = clinical_data(:,8);
site = clinical_data(:,2);

covariates = [ptsd_dx age sex ptsd_severity group index_trauma];
tic
result_step_two = stack_sub_matrices_function_combat_enigma(sub_ids, good_ROIs, resting_filename, ...,
    atlas_location,mask_values,site,covariates);
save result_step_two result_step_two
toc

big_r=result_step_two.big_r;
all_sub_r = result_step_two.all_sub_r;
big_r_resid = result_step_two.big_r_resid;
all_sub_r_resid = result_step_two.harmonized;
    

%% next step: can only have positive connections, so take negatives out
pos_big_r=big_r_resid;
pos_big_r(pos_big_r<0)=0;

%need to have sparse network, so use this function to identify the maximum
%threshold that results in a fully connected network
thresh_z=calculate_max_connected_threshold(pos_big_r);

%now threshold the positive network so that only connections above this
%threshold are retained
pos_big_r_thresh=pos_big_r; 
pos_big_r_thresh(pos_big_r_thresh<thresh_z)=0;
%calculate density of the network just for characterization
d=density_und(pos_big_r_thresh);

num_iter=10000;
spatial_map_filename=['name_of_group_level_map_' num2str(num_iter) '.nii'];
tic
[message result_step_four] = calculate_group_level_community_structure_function_surprise_500(atlas_location, num_iter, pos_big_r_thresh, spatial_map_filename,mask_values);
toc
save result_step_four result_step_four

community_structure=result_step_four.community_structure;
community_sort_index=result_step_four.community_sort_index;
communities = numel(unique(community_structure));


%% Call function to do large-scale network within-subject analysis 
% determines strength of coupling between networks within an individual

subject_spatial_map_output_prefix = ['/path/to/and/name/of/subject/maps'];

%make sa for loop that creates a string array for sub_ids, sub_num
%to be user-defined

clear sub_num
for sn=1:numel(sub_ids)
    
    if sub_ids(sn)<=41 %Tours
        if sub_ids(sn)<10
            sub_num{sn}=['000' num2str(sub_ids(sn))];
        elseif sub_ids(sn)>=10&&sub_ids(sn)<100
            sub_num{sn}=['00' num2str(sub_ids(sn))];
        end
    
    elseif sub_ids(sn) >= 42 && sub_ids(sn) <= 126 %Emory
        if sub_ids(sn)>=10&&sub_ids(sn)<100
            sub_num{sn}=['00' num2str(sub_ids(sn))];
        elseif sub_ids(sn) >=100 
            sub_num{sn}=['0' num2str(sub_ids(sn))];
        end
    
    elseif sub_ids(sn)  >= 127 && sub_ids(sn) <=206 %McLean 
        sub_num{sn}=['0' num2str(sub_ids(sn))];
   
    elseif sub_ids(sn) >= 207 && sub_ids(sn) <=256 %Vanderbilt
        sub_num{sn}=['0' num2str(sub_ids(sn))];
    
    elseif sub_ids(sn) >= 257 && sub_ids(sn) <=362 %UWCisler
        sub_num{sn}=['0' num2str(sub_ids(sn))];
 
    elseif sub_ids(sn) >= 363 && sub_ids(sn) <=395 %Gronigen
        sub_num{sn}=['0' num2str(sub_ids(sn))];    

    elseif sub_ids(sn) >= 445 && sub_ids(sn) <=507 %Michigan
        sub_num{sn}=['0' num2str(sub_ids(sn))]; 
                
    elseif sub_ids(sn) >= 508 && sub_ids(sn) <=551 %Munster
        sub_num{sn}=['0' num2str(sub_ids(sn))];     
    
    elseif sub_ids(sn) >= 552 && sub_ids(sn) <=633 %Duke
        sub_num{sn}=['0' num2str(sub_ids(sn))];      
    
    elseif sub_ids(sn) >= 634 && sub_ids(sn) <=673 %UMN
        sub_num{sn}=['0' num2str(sub_ids(sn))];       

    elseif sub_ids(sn) >= 674 && sub_ids(sn) <=731 %UWGrupe
        sub_num{sn}=['0' num2str(sub_ids(sn))];      
    
    elseif sub_ids(sn) >= 732 && sub_ids(sn) <=790 %AMC
        sub_num{sn}=['0' num2str(sub_ids(sn))];     
    
    elseif sub_ids(sn) >= 791 && sub_ids(sn) <=868 %Milwaukee
        sub_num{sn}=['0' num2str(sub_ids(sn))];
    
    elseif sub_ids(sn) >= 869 && sub_ids(sn) <=927 %Utrecht
        sub_num{sn}=['0' num2str(sub_ids(sn))];
    
    elseif sub_ids(sn) >= 928 && sub_ids(sn) <=999 %Minneapolis VA
        sub_num{sn}=['0' num2str(sub_ids(sn))];
       
    elseif sub_ids(sn) >= 1000 && sub_ids(sn) <=1070 %Stanford Brains
        sub_num{sn}=num2str(sub_ids(sn)); 
        
    elseif sub_ids(sn) >= 1071 && sub_ids(sn) <=1211 %Western Ontario
        sub_num{sn}=num2str(sub_ids(sn)); 
    
    elseif sub_ids(sn) >= 1212 && sub_ids(sn) <=1271  %Stanford CausCon
        sub_num{sn}=num2str(sub_ids(sn));
    
    elseif sub_ids(sn) >= 1272 && sub_ids(sn) <=1334  %Masaryk
        sub_num{sn}=num2str(sub_ids(sn));
    
    elseif sub_ids(sn) >= 1335 && sub_ids(sn) <=1408  %Columbia
        sub_num{sn}=num2str(sub_ids(sn));
   
    elseif sub_ids(sn) >= 1409 && sub_ids(sn) <=1459  %Ghent
        sub_num{sn}=num2str(sub_ids(sn));
   
    elseif sub_ids(sn) >= 1460 && sub_ids(sn) <= 1499  %Toledo MVA
        sub_num{sn}=num2str(sub_ids(sn));
    
    elseif sub_ids(sn) >= 1500 && sub_ids(sn) <=1526  %MToledo ONG
        sub_num{sn}=num2str(sub_ids(sn));
    
    elseif sub_ids(sn) >= 1537 %Waco
        sub_num{sn}=num2str(sub_ids(sn)); 
    end
end
% do thresholding for all_sub_r

pos_all_sub_r=all_sub_r_resid;
pos_all_sub_r(pos_all_sub_r<0)=0;

thresh_z=zeros(size(pos_all_sub_r,3),1);
pos_all_sub_r_thresh=zeros(size(pos_all_sub_r));
tic
for sn=1:size(pos_all_sub_r,3)
    thresh_z(sn,1)=calculate_max_connected_threshold(squeeze(pos_all_sub_r(:,:,sn)));
    temp_r=squeeze(pos_all_sub_r(:,:,sn));
    temp_r(temp_r<thresh_z(sn,1))=0;
    pos_all_sub_r_thresh(:,:,sn)=temp_r;
end
toc

save pos_all_sub_r_thresh250 pos_all_sub_r_thresh

%call function to apply PACO-defined group-level map to individual subjects

pre_defined_network_location='path/to/group-level/spatial_map.nii';
%rename pos_all_sub_r_thresh to all_sub_r so the function will take it
all_sub_r=pos_all_sub_r_thresh;

tic
result_step_five_surprise_sample_mask = calculate_pre_defined_network_within_subject(atlas_location,sub_num,...,
    pre_defined_network_location, subject_spatial_map_output_prefix,community_structure,mask_values,all_sub_r);
toc
save result_step_five_surprise_sample_mask_250_full result_step_five_surprise_sample_mask



