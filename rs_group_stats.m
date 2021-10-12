%%
%Author: Marisa C. Ross, PhD
%Date: May 2021

addpath(genpath('/Volumes/Vol2/cisler/DOP/network/'))
addpath(genpath('/Volumes/Vol2/cisler/matlab_toolboxes/'))
addpath(genpath('/Volumes/Vol2/cisler/network/functions/'))
addpath(genpath('/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/final_analyses/BrainNetViewer_20191031/'))
clear

%%
%define initial variables  

sub_ids = load('/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/site_testing/sub_ids_all_col2');


%% Load in and Define Clinical Data
clinical_data=load(['/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/site_testing/clinical_data_all_col2']);
keep_clinical=find(ismember(clinical_data(:,1),sub_ids)==1);
clinical_data=clinical_data(keep_clinical,:);

age = clinical_data(:,3); 
ptsd_severity=clinical_data(:,7); 
sex = clinical_data(:,4); %female = -.5, male = 0.5
ptsd_dx = clinical_data(:,5); %control = -0.5, ptsd = 0.5
group = clinical_data(:,6); %-1 = NTC, 0 = TEC, 1 = PTSD
index_trauma = clinical_data(:,8);
site = clinical_data(:,2);
subs = clinical_data(:,1);

%% Load in network results
load result_step_one_all_col2.mat
good_ROIs = result_step_one_surprise_250.good_ROIs;
mask_values = result_step_one_surprise_250.mask_values;

load result_step_two_all_col2.mat
group_r = result_step_two.big_r_resid;
all_sub_r = result_step_two.harmonized;

load result_step_four_all_col2.mat
community_structure = result_step_four.community_structure;
community_sort_index = result_step_four.community_sort_index;
communities_group = numel(unique(community_structure));
surprise = result_step_four.surprise;

load result_step_five_full_sample.mat
pos_par = result_step_five_surprise_sample_mask.all_graph_theory_indices(:,:,4);
roi_strength = result_step_five_surprise_sample_mask.all_graph_theory_indices(:,:,2);
mean_pos_par = mean(pos_par,2);
mean_roi_strength = mean(roi_strength,2);

load all_whole_brain_strength 
load all_cen_strength
load all_dmn_strength
load mean_pos_par_dmn 
load mean_pos_par_cen
load minus_dmn_pos_par
load minus_dmn_strength
load minus_cen_pos_par
load minus_cen_strength


%% Question 1: does network org differ by group, and does binary sex moderate the relationship between group and network org?

%create strings for dummy coding group
group_str=string(group);
group_str(group_str=="-1")="NTC";
group_str(group_str=="0")="TEC";
group_str(group_str=="1")="PTSD";
group_cat = categorical(group_str);
group_dummy = dummyvar(group_cat);
group_dummy(:,2) = []; %deletes the PTSD column to create a reference group; dv1 is NTC vs. PTSD, dv2 is TEC vs. PTSD

%mean participation coefficient
Y=mean_pos_par;
X=[group_dummy sex age];
random_X=[site];
stack_data=[Y X random_X subs];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','group1','group2','sex','age','site','sub'});
lme1=fitlme(tbl,'mean_pos_par~group1*sex+group2*sex+age+(1+site|sub)','DummyVarCoding','reference');

%DMN
Y=mean_pos_par_dmn;
X=[group_dummy sex age minus_dmn_pos_par];
random_X=[site];
stack_data=[Y X random_X subs];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','group1','group2','sex','age','wb','site','sub'});
lme1_1=fitlme(tbl,'mean_pos_par~group1*sex+group2*sex+age+wb+(1+site|sub)','DummyVarCoding','reference');

%CEN
Y=mean_pos_par_cen;
X=[group_dummy sex age minus_cen_pos_par];
random_X=[site];
stack_data=[Y X random_X subs];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','group1','group2','sex','age','wb','site','sub'});
lme1_2=fitlme(tbl,'mean_pos_par~group1*sex+group2*sex+age+wb+(1+site|sub)','DummyVarCoding','reference');

%mean connectivity strength across the brain
Y=all_whole_brain_strength;
X=[group_dummy sex age];
random_X=[site];
stack_data=[Y X random_X subs];
tbl=array2table(stack_data,'variablenames',{'strength','group1','group2','sex','age','site','sub'});
lme1b=fitlme(tbl,'strength~group1*sex+group2*sex+age+(1+site|sub)','DummyVarCoding','reference');

%DMN
Y=all_dmn_strength;
X=[group_dummy sex age minus_dmn_strength];
random_X=[site];
stack_data=[Y X random_X subs];
tbl=array2table(stack_data,'variablenames',{'strength','group1','group2','sex','age','wb','site','sub'});
lme1b_1=fitlme(tbl,'strength~group1*sex+group2*sex+age+wb+(1+site|sub)','DummyVarCoding','reference');

%CEN
Y=all_cen_strength;
X=[group_dummy sex age minus_cen_strength];
random_X=[site];
stack_data=[Y X random_X subs];
tbl=array2table(stack_data,'variablenames',{'strength','group1','group2','sex','age','wb','site','sub'});
lme1b_2=fitlme(tbl,'strength~group1*sex+group2*sex+age+wb+(1+site|sub)','DummyVarCoding','reference');

%make a figure for the main effect of group in pos_par

%get residuals
stats = regstats(mean_pos_par,[sex age]);stats.tstat.t
pos_par_resid = mean_pos_par+stats.standres;

group_means(:,1) = mean(pos_par_resid(group==-1));
group_means(:,2) = mean(pos_par_resid(group==0));
group_means(:,3) = mean(pos_par_resid(group==1));

group_se(:,1) = std(pos_par_resid(group==-1))./sqrt(numel(group==-1));
group_se(:,2) = std(pos_par_resid(group==0))./sqrt(numel(group==0));
group_se(:,3) = std(pos_par_resid(group==1))./sqrt(numel(group==1));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [.627 .627 .627];
hbar.CData(2,:) = [.4 .689 1];
hbar.CData(3,:) = [1 .4 .4];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Whole-Brain Average Participation Coefficient');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('group_WB_pos_par','-dtiff')

%make a figure for the main effect of sex in pos_par
stats = regstats(mean_pos_par,[group age]);stats.tstat.t
pos_par_resid = mean_pos_par+stats.standres;

group_means(:,1) = mean(pos_par_resid(sex==-0.5));
group_means(:,2) = mean(pos_par_resid(sex==0.5));

group_se(:,1) = std(pos_par_resid(sex==-0.5))./sqrt(numel(sex==-0.5));
group_se(:,2) = std(pos_par_resid(sex==0.5))./sqrt(numel(sex==0.5));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [.60 0.2 1];
hbar.CData(2,:) = [.88 .88 .88];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Whole-Brain Average Participation Coefficient');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('ME_sex_pos_par','-dtiff')

%make a figure for the main effect of group in dmn_pos_par

%get residuals
stats = regstats(mean_pos_par_dmn,[sex age minus_dmn_pos_par]);stats.tstat.t
pos_par_resid = mean_pos_par_dmn+stats.standres;

group_means(:,1) = mean(pos_par_resid(group==-1));
group_means(:,2) = mean(pos_par_resid(group==0));
group_means(:,3) = mean(pos_par_resid(group==1));

group_se(:,1) = std(pos_par_resid(group==-1))./sqrt(numel(group==-1));
group_se(:,2) = std(pos_par_resid(group==0))./sqrt(numel(group==0));
group_se(:,3) = std(pos_par_resid(group==1))./sqrt(numel(group==1));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [.627 .627 .627];
hbar.CData(2,:) = [.4 .689 1];
hbar.CData(3,:) = [1 .4 .4];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Average DMN Participation Coefficient');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('group_DMN_pos_par','-dtiff')

%make a figure for the main effect of group in cen_pos_par

%get residuals
stats = regstats(mean_pos_par_cen,[sex age minus_cen_pos_par]);stats.tstat.t
pos_par_resid = mean_pos_par_cen+stats.standres;

group_means(:,1) = mean(pos_par_resid(group==-1));
group_means(:,2) = mean(pos_par_resid(group==0));
group_means(:,3) = mean(pos_par_resid(group==1));

group_se(:,1) = std(pos_par_resid(group==-1))./sqrt(numel(group==-1));
group_se(:,2) = std(pos_par_resid(group==0))./sqrt(numel(group==0));
group_se(:,3) = std(pos_par_resid(group==1))./sqrt(numel(group==1));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [.627 .627 .627];
hbar.CData(2,:) = [.4 .689 1];
hbar.CData(3,:) = [1 .4 .4];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Average CEN Participation Coefficient');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('group_CEN_pos_par','-dtiff')

%make a figure for the main effect of group in wb_strength

%get residuals
stats = regstats(all_whole_brain_strength,[sex age]);stats.tstat.t
strength_resid = all_whole_brain_strength+stats.standres;

group_means(:,1) = mean(strength_resid(group==-1));
group_means(:,2) = mean(strength_resid(group==0));
group_means(:,3) = mean(strength_resid(group==1));

group_se(:,1) = std(strength_resid(group==-1))./sqrt(numel(group==-1));
group_se(:,2) = std(strength_resid(group==0))./sqrt(numel(group==0));
group_se(:,3) = std(strength_resid(group==1))./sqrt(numel(group==1));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
%set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [.627 .627 .627];
hbar.CData(2,:) = [.4 .689 1];
hbar.CData(3,:) = [1 .4 .4];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Whole-Brain Average Connectivity Strength');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('group_WB_strength','-dtiff')

%make a figure for the main effect of group in dmn_strength

%get residuals
stats = regstats(all_dmn_strength,[sex age minus_dmn_strength]);stats.tstat.t
strength_resid = all_dmn_strength+stats.standres;

group_means(:,1) = mean(strength_resid(group==-1));
group_means(:,2) = mean(strength_resid(group==0));
group_means(:,3) = mean(strength_resid(group==1));

group_se(:,1) = std(strength_resid(group==-1))./sqrt(numel(group==-1));
group_se(:,2) = std(strength_resid(group==0))./sqrt(numel(group==0));
group_se(:,3) = std(strength_resid(group==1))./sqrt(numel(group==1));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
%set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [.627 .627 .627];
hbar.CData(2,:) = [.4 .689 1];
hbar.CData(3,:) = [1 .4 .4];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Within DMN Connectivity Strength');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('group_DMN_strength','-dtiff')

%make a figure for the main effect of group in cen_pos_par

%get residuals
stats = regstats(all_cen_strength,[sex age minus_cen_strength]);stats.tstat.t
strength_resid = all_cen_strength+stats.standres;

group_means(:,1) = mean(strength_resid(group==-1));
group_means(:,2) = mean(strength_resid(group==0));
group_means(:,3) = mean(strength_resid(group==1));

group_se(:,1) = std(strength_resid(group==-1))./sqrt(numel(group==-1));
group_se(:,2) = std(strength_resid(group==0))./sqrt(numel(group==0));
group_se(:,3) = std(strength_resid(group==1))./sqrt(numel(group==1));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
%set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [.627 .627 .627];
hbar.CData(2,:) = [.4 .689 1];
hbar.CData(3,:) = [1 .4 .4];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Within CEN Connectivity Strength');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('group_CEN_strength','-dtiff')

%make a figure for the main effect of sex in dmn strength
stats = regstats(all_dmn_strength,[group age minus_dmn_strength]);stats.tstat.t
pos_par_resid = all_dmn_strength+stats.standres;

group_means(:,1) = mean(pos_par_resid(sex==-0.5));
group_means(:,2) = mean(pos_par_resid(sex==0.5));

group_se(:,1) = std(pos_par_resid(sex==-0.5))./sqrt(numel(sex==-0.5));
group_se(:,2) = std(pos_par_resid(sex==0.5))./sqrt(numel(sex==0.5));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [.60 0.2 1];
hbar.CData(2,:) = [.88 .88 .88];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Within DMN Connectivity Strength');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('ME_sex_dmn_strength','-dtiff')

%% Question 2: Does network org differ with symptom severity in trauma exposed groups? Does binary sex moderate this relationship?

%mean participation coefficient across the brain
Y=mean_pos_par(group==0 | group==1);
X=[ptsd_severity(group==0 | group==1) sex(group==0 | group==1) age(group==0 | group==1)];
random_X=[site(group==0 | group==1)];
stack_data=[Y X random_X subs(group==0 | group==1)];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','sev','sex','age','site','sub'});
lme2=fitlme(tbl,'mean_pos_par~sev*sex+age+(1+site|sub)');

%mean participation coefficient of DMN
Y=mean_pos_par_dmn(group==0 | group==1);
X=[ptsd_severity(group==0 | group==1) sex(group==0 | group==1) age(group==0 | group==1) minus_dmn_pos_par(group==0 | group==1)];
random_X=[site(group==0 | group==1)];
stack_data=[Y X random_X subs(group==0 | group==1)];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','sev','sex','age','wb','site','sub'});
lme2_1=fitlme(tbl,'mean_pos_par~sev*sex+age+wb+(1+site|sub)');

%mean participation coefficient of CEN
Y=mean_pos_par_cen(group==0 | group==1);
X=[ptsd_severity(group==0 | group==1) sex(group==0 | group==1) age(group==0 | group==1) minus_cen_pos_par(group==0 | group==1)];
random_X=[site(group==0 | group==1)];
stack_data=[Y X random_X subs(group==0 | group==1)];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','sev','sex','age','wb','site','sub'});
lme2_2=fitlme(tbl,'mean_pos_par~sev*sex+age+wb+(1+site|sub)');

%mean connectivity strength across the brain
Y=all_whole_brain_strength(group==0 | group==1);
X=[ptsd_severity(group==0 | group==1) sex(group==0 | group==1) age(group==0 | group==1)];
random_X=[site(group==0 | group==1)];
stack_data=[Y X random_X subs(group==0 | group==1)];
tbl=array2table(stack_data,'variablenames',{'strength','sev','sex','age','site','sub'});
lme2b=fitlme(tbl,'strength~sev*sex+age+(1+site|sub)');

%mean connectivity strength DMN
Y=all_dmn_strength(group==0 | group==1);
X=[ptsd_severity(group==0 | group==1) sex(group==0 | group==1) age(group==0 | group==1) minus_dmn_strength(group==0 | group==1)];
random_X=[site(group==0 | group==1)];
stack_data=[Y X random_X subs(group==0 | group==1)];
tbl=array2table(stack_data,'variablenames',{'strength','sev','sex','age','wb','site','sub'});
lme2b_1=fitlme(tbl,'strength~sev*sex+age+wb+(1+site|sub)');

%mean connectivity strength CEN
Y=all_cen_strength(group==0 | group==1);
X=[ptsd_severity(group==0 | group==1) sex(group==0 | group==1) age(group==0 | group==1) minus_cen_strength(group==0 | group==1)];
random_X=[site(group==0 | group==1)];
stack_data=[Y X random_X subs(group==0 | group==1)];
tbl=array2table(stack_data,'variablenames',{'strength','sev','sex','age','wb','site','sub'});
lme2b_1=fitlme(tbl,'strength~sev*sex+age+wb+(1+site|sub)');


%% Question 3: Does index trauma type moderate the relationship between ptsd severity and network org?

%create strings for dummy coding index trauma
index_str=string(index_trauma(index_trauma < 99));
index_str(index_str=="-1")="COM";
index_str(index_str=="0")="IPV";
index_str(index_str=="1")="OT";
index_cat = categorical(index_str);
index_dummy = dummyvar(index_cat);
index_dummy(:,2) = []; %deletes the IPV column to create a reference group; dv1 is COM vs. IPV, dv2 is OT vs. IPV

%average  participation coefficient
Y=mean_pos_par(index_trauma<99);
X=[ptsd_severity(index_trauma<99) index_dummy sex(index_trauma<99) age(index_trauma<99)];
random_X=site(index_trauma<99);
stack_data=[Y X random_X subs(index_trauma<99)];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','sev','dv1','dv2','sex','age','site','sub'});
lme3=fitlme(tbl,'mean_pos_par~sev*dv1+sev*dv2+sex+age+(1+site|sub)','DummyVarCoding','reference');

%average  participation coefficient DMN
Y=mean_pos_par_dmn(index_trauma<99);
X=[ptsd_severity(index_trauma<99) index_dummy sex(index_trauma<99) age(index_trauma<99) minus_dmn_pos_par(index_trauma<99)];
random_X=site(index_trauma<99);
stack_data=[Y X random_X subs(index_trauma<99)];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','sev','dv1','dv2','sex','age','wb','site','sub'});
lme3_1=fitlme(tbl,'mean_pos_par~sev*dv1+sev*dv2+sex+age+wb+(1+site|sub)','DummyVarCoding','reference');

%average  participation coefficient CEN
Y=mean_pos_par_cen(index_trauma<99);
X=[ptsd_severity(index_trauma<99) index_dummy sex(index_trauma<99) age(index_trauma<99) minus_cen_pos_par(index_trauma<99)];
random_X=site(index_trauma<99);
stack_data=[Y X random_X subs(index_trauma<99)];
tbl=array2table(stack_data,'variablenames',{'mean_pos_par','sev','dv1','dv2','sex','age','wb','site','sub'});
lme3_2=fitlme(tbl,'mean_pos_par~sev*dv1+sev*dv2+sex+age+wb+(1+site|sub)','DummyVarCoding','reference');

%average whole-brain strength
Y=all_whole_brain_strength(index_trauma<99);
X=[ptsd_severity(index_trauma<99) index_dummy sex(index_trauma<99) age(index_trauma<99)];
random_X=site(index_trauma<99);
stack_data=[Y X random_X subs(index_trauma<99)];
tbl=array2table(stack_data,'variablenames',{'strength','sev','dv1','dv2','sex','age','site','sub'});
lme3b=fitlme(tbl,'strength~sev*dv1+sev*dv2+sex+age+(1+site|sub)','DummyVarCoding','reference');

%average DMN strength
Y=all_dmn_strength(index_trauma<99);
X=[ptsd_severity(index_trauma<99) index_dummy sex(index_trauma<99) age(index_trauma<99) minus_dmn_strength(index_trauma<99)];
random_X=site(index_trauma<99);
stack_data=[Y X random_X subs(index_trauma<99)];
tbl=array2table(stack_data,'variablenames',{'strength','sev','dv1','dv2','sex','age','wb','site','sub'});
lme3b_1=fitlme(tbl,'strength~sev*dv1+sev*dv2+sex+age+wb+(1+site|sub)','DummyVarCoding','reference');

%average CEN strength
Y=all_cen_strength(index_trauma<99);
X=[ptsd_severity(index_trauma<99) index_dummy sex(index_trauma<99) age(index_trauma<99) minus_cen_strength(index_trauma<99)];
random_X=site(index_trauma<99);
stack_data=[Y X random_X subs(index_trauma<99)];
tbl=array2table(stack_data,'variablenames',{'strength','sev','dv1','dv2','sex','age','wb','site','sub'});
lme3b_2=fitlme(tbl,'strength~sev*dv1+sev*dv2+sex+age+wb+(1+site|sub)','DummyVarCoding','reference');

%graph the group effect in CEN pos par
stats = regstats(mean_pos_par_cen,[sex age minus_cen_pos_par]);stats.tstat.t
pos_par_resid = mean_pos_par_cen+stats.standres;

group_means(:,1) = mean(pos_par_resid(index_trauma==-1));
group_means(:,2) = mean(pos_par_resid(index_trauma==0));
group_means(:,3) = mean(pos_par_resid(index_trauma==1));

group_se(:,1) = std(pos_par_resid(index_trauma==-1))./sqrt(numel(index_trauma==-1));
group_se(:,2) = std(pos_par_resid(index_trauma==0))./sqrt(numel(index_trauma==0));
group_se(:,3) = std(pos_par_resid(index_trauma==1))./sqrt(numel(index_trauma==1));

hbar=bar([group_means'],1);
hbar.FaceColor='flat';
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
set(hbar,'EdgeColor','k','LineWidth',1.5)
%set(gca,'ytick',-1:.2:1, 'yticklabels',-1:.2:1)
hbar.CData(1,:)  = [0 .6 .298];
hbar.CData(2,:) = [.6 .6 0];
hbar.CData(3,:) = [.0 .6 .6];
set(gca,'FontWeight','bold','Fontsize',10)
set(gca,'xticklabels',[]);
this_y_label=ylabel('Residualized Within CEN Participation Coefficient');
hold on
errorbar(5:1:1.5,group_means(1,:),group_se(1,:),'.k','MarkerSize',1,'LineWidth',1.5)
print('index_trauma_CEN_strength','-dtiff')


%% ROI-by-ROI break down of graph-theory indices

% First, make .nii of good_ROIs so number of nodes is accurate

ROIs = load_nii('/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/full_analyses/all_subs_good_ROIs_222.nii');
ROI_size = size(ROIs.img);
in_brain = find(ROIs.img > 0);
%in_brain = unique(ROIs.img(in_brain));

store_all_pos_par=zeros(numel(sub_ids),numel(in_brain));

for sn=1:numel(sub_ids)
    if sub_ids(sn)<10
       sub_num=['000' num2str(sub_ids(sn))];
    elseif sub_ids(sn)>=10&&sub_ids(sn)<100
       sub_num=['00' num2str(sub_ids(sn))];
    elseif sub_ids(sn)>=100&&sub_ids(sn)<1000
       sub_num=['0' num2str(sub_ids(sn))];
    elseif sub_ids(sn)>1000
        sub_num=num2str(sub_ids(sn));
    end
    disp(sub_num)
    data=load_nii(['/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/full_analyses/sub_maps/rest.surprise_250_' sub_num '.nii']);
    sub_brain=double(reshape(squeeze(data.img),ROI_size(1)*ROI_size(2)*ROI_size(3),size(squeeze(data.img),4)));
    select_briks = size(sub_brain,3):size(sub_brain,2);
    sub_briks = sub_brain(in_brain,select_briks);
    store_all_briks(sn,:,:,:) = sub_briks;
    sub_brain_pos_par=sub_brain(in_brain,4); % this saves the brik containing the pos_par for every subject, for every node, in a vector
    store_all_pos_par(sn,:)=sub_brain_pos_par;
end

%to use whole_brain_regression function
%make design matrix
X=[ptsd_dx sex age site];
filename='ptsd_v_control_t_ROI_250_resid';
whole_brain_data=store_all_pos_par;%data you care about in participant x brain space
mask=in_brain;
%data is an original nifti file loaded in and untouched
message = whole_brain_regression(whole_brain_data,X,filename,mask,data);

% Node analysis at the group level to determine critical t-values
% load in all relevant inputs

all_graph_theory_indices_resid=zeros(numel(sub_ids),numel(mask_values),7);
for l=1:numel(mask_values)
    find_roi=find(ROIs.img(in_brain)==mask_values(l));
    temp=squeeze(store_all_briks(:,find_roi,:));
    all_graph_theory_indices_resid(:,l,:) = squeeze(temp(:,1,:));
end
save all_graph_theory_indices all_graph_theory_indices_resid

atlas_location = '/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/full_analyses/all_subs_good_ROIs_222.nii';
covariates = [age sex site];
clinical_variables = ptsd_dx;
clinical_labels = {'ptsd'};
group_spatial_map_output_prefix = '/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/full_analyses/ptsd_v_control_map_';
p_crit = 0.05;
stat_choice = 2; %stat_choice 1 = partial correlation (spearman), stat_choice 2 = regression (tstat)

result = group_level_node_analysis(all_graph_theory_indices_resid,atlas_location,covariates,clinical_variables,clinical_labels,mask_values,group_spatial_map_output_prefix,p_crit,stat_choice);

save group_level_node_result result


%% Node-level group analysis

load group_level_node_result.mat
critcal_values = stat_crit_values;
t_values = result.save_stat;
p_values = result.save_ps;

% order of graph theory indices and briks in ptsd_v_control_ptsd.nii with critical t values:
% [community_structure (1.0396) roi_strength (3.1377) roi_clust (3.3042) pos_par (3.3919) neg_par(nan) ...,
% roi_within_z (nan) roi_betweeness (3.5029)]

