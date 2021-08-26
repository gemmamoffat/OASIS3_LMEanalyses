%% categorize participants into progressing vs non-progressing categories
clindata=readtable('clinical_longitudinal.csv'); 
clindata= sortrows(clindata,[2 3]);
allsubs=unique(clindata.Subject);

for i=1:length(unique(clindata.Subject))
    id=strcmp(clindata.Subject, allsubs{i}); 
    cdr=clindata.cdr(id); mmse=clindata.mmse(id); date=clindata.DaysSinceBaseline(id)/365;
    X=[ones(length(cdr),1), date];
    timeline(i)=max(clindata.DaysSinceBaseline(id))/365;
    if length(cdr) > 1
        b=regress(cdr,X); c=regress(mmse,X);
        cdrslopes(i)=b(2); mmseslopes(i)=c(2);
        if b(2)>0.0001
             decrease=0;
             for j=1:length(cdr)-2
                 if cdr(j) > cdr(j+1) 
                     if length(cdr) > 2
                         if cdr(j+2)==cdr(j+1)
                            decrease=1;
                         end
                     else
                         decrease=1;
                     end
                 end
             end
             if ((decrease == 0) || cdr(length(cdr))-cdr(1)>1)
                delta(i)=1;
             else 
                 if (mmseslopes(i) < 0)
                     delta(i)=1;
                 else
                    delta(i)=0;
                 end 
             end
        else
            delta(i)=0;
        end
    else
        cdrslopes(i)=NaN;
        mmseslopes(i)=NaN;
    end
    dx1(i)=sum(clindata.dx1_num(id));
    apoe(i)=mean(clindata.apoe(id));
    baseline_cdr(i)=cdr(1);
    age(i)=mean(clindata.ageAtEntry(id));
    sex(i)=clindata.M_F{id,1};
    education(i)=mean(clindata.Education(id));
    race(i)=unique(clindata.Race(id,1));
    years_in_study(i)=max(clindata.DaysSinceBaseline(id)/365);
end

%% import amyloid data; add demographic variables for all amyloid participants
amyloid=readtable('amyloid_longitudinal.csv');
for i=1:length(amyloid.Subject)  
    if ismember(amyloid.Subject{i},allsubs)
        id=strcmp(amyloid.Subject{i}, allsubs); 
        d=delta(id==1); 
        amyloid.delta(i,1)=d;clear d; 
        amyloid.baseline_age(i,1)=age(id==1);
        amyloid.baseline_cdr(i,1)=baseline_cdr(id==1);
        amyloid.age(i,1)=amyloid.time_in_y(i)+age(id==1);
        amyloid.sex(i,1)=sex(id==1);
        amyloid.education(i,1)=education(id==1);
        amyloid.race(i,1)=race(id==1);
        amyloid.apoe(i,1)=apoe(id==1);
        amyloid.dx1(i,1)=dx1(id==1);
    end   
end

allsubs_staticdata=clindata(:,{'Subject','ageAtEntry','Education','Race','Ethnicity','M_F','cdr'});    
for i=1:length(clindata.Subject)
    id=strcmp(clindata.Subject{i}, allsubs);
    allsubs_staticdata.delta(i,1)=delta(id==1);
    allsubs_staticdata.cdrslope(i,1)=cdrslopes(id==1);
    allsubs_staticdata.baseline_cdr(i,1)=baseline_cdr(id==1);
    allsubs_staticdata.years_in_study(i,1)=years_in_study(id==1);
    allsubs_staticdata.mmseslope(i,1)=mmseslopes(id==1);
    allsubs_staticdata.apoe(i,1)=apoe(id==1);
end

%% exculsions

agem_amyloid=amyloid(:,:);
agem_amyloid(agem_amyloid.baseline_age<57 & agem_amyloid.delta==0,:)=[];
agem_amyloid(agem_amyloid.dx1>0,:)=[];
agem_repeated_amyloid=grouptransform(agem_amyloid,'Subject',@numel,'ReplaceValues', false); % counts how many measures per subject, adds column of vaues
agem_repeated_amyloid=agem_repeated_amyloid(1:length(agem_repeated_amyloid.Subject),1:286); % only keeps valid columns
agem_repeated_amyloid.Properties.VariableNames{286} = 'measures';
agem_repeated_amyloid(agem_repeated_amyloid.measures==1,:)=[];
agem_repeated_amyloid = sortrows(agem_repeated_amyloid,[2 6]); % orders rows based on time_in_y after subject ID


% make table with one row per participant
[C,ia,ic] = unique(agem_repeated_amyloid.Subject); 
agem_unique_repeated_amyloid = agem_repeated_amyloid(ia,:); 

clear C; clear ia; clear ic;
[C,ia,ic] = unique(allsubs_staticdata.Subject); 
agem_subject_stats = allsubs_staticdata(ia,:);

agem_subject_stats(~ismember(agem_subject_stats.Subject,agem_repeated_amyloid.Subject),:)=[];

for i=1:length(agem_unique_repeated_amyloid.Subject)
    if agem_unique_repeated_amyloid.time_in_y(i) > agem_subject_stats.years_in_study(i)
        agem_subject_stats.years_in_study(i)=agem_unique_repeated_amyloid.time_in_y(i);
    end
end

% figure;
% histogram(agem_unique_repeated_amyloid.delta,'FaceColor',[0.6 0.7 0.9],'EdgeColor',[0.6 0.8 1]);
% title('Subject classification based on change in cognitive status over time','FontSize',15);
% xlabel('delta','FontSize',15);
% ylabel('Number of subjects','FontSize',15);

figure;
v=violinplot(agem_subject_stats.cdrslope, agem_subject_stats.delta, 'ShowData', logical(0), 'ShowMean', logical(1),'BoxColor',[0.6 0.7 0.9]);
figure;
v=violinplot(agem_subject_stats.mmseslope, agem_subject_stats.delta, 'ShowData', logical(0), 'ShowMean', logical(1),'BoxColor',[0.6 0.7 0.9]);

agem_subject_stats=agem_subject_stats(ismember(agem_subject_stats.Subject,agem_unique_repeated_amyloid.Subject),:);
sortrows(agem_subject_stats,[8 1]);

%% make table 1

table1{1,1}='Subjects';
table1{1,2}=sum(agem_subject_stats.delta==0);
table1{1,3}=sum(agem_subject_stats.delta==1);
table1{2,1}='Age';
table1{2,2}=mean(agem_subject_stats.ageAtEntry(agem_subject_stats.delta==0));
table1{2,3}=mean(agem_subject_stats.ageAtEntry(agem_subject_stats.delta==1));
table1{2,4}=std(agem_subject_stats.ageAtEntry(agem_subject_stats.delta==0));
table1{2,5}=std(agem_subject_stats.ageAtEntry(agem_subject_stats.delta==1));
table1{2,6}=anova1(agem_subject_stats.ageAtEntry,agem_subject_stats.delta);
table1{3,1}='Gender - %Male';
table1{3,2}=sum(strcmp(agem_subject_stats.M_F(agem_subject_stats.delta==0),'M'));
table1{3,3}=sum(strcmp(agem_subject_stats.M_F(agem_subject_stats.delta==1),'M'));
table1{4,1}='Education';
table1{4,2}=mean(agem_subject_stats.Education(agem_subject_stats.delta==0));
table1{4,3}=mean(agem_subject_stats.Education(agem_subject_stats.delta==1));
table1{4,4}=std(agem_subject_stats.Education(agem_subject_stats.delta==0));
table1{4,5}=std(agem_subject_stats.Education(agem_subject_stats.delta==1));
table1{4,6}=anova1(agem_subject_stats.Education,agem_subject_stats.delta);
table1{5,1}='CDR slope';
table1{5,2}=mean(agem_subject_stats.cdrslope(agem_subject_stats.delta==0));
table1{5,3}=mean(agem_subject_stats.cdrslope(agem_subject_stats.delta==1));
table1{5,4}=std(agem_subject_stats.cdrslope(agem_subject_stats.delta==0));
table1{5,5}=std(agem_subject_stats.cdrslope(agem_subject_stats.delta==1));
table1{5,6}=anova1(agem_subject_stats.cdrslope,agem_subject_stats.delta);
table1{6,1}='MMSE slope';
table1{6,2}=mean(agem_subject_stats.mmseslope(agem_subject_stats.delta==0));
table1{6,3}=mean(agem_subject_stats.mmseslope(agem_subject_stats.delta==1));
table1{6,4}=std(agem_subject_stats.mmseslope(agem_subject_stats.delta==0));
table1{6,5}=std(agem_subject_stats.mmseslope(agem_subject_stats.delta==1));
table1{6,6}=anova1(agem_subject_stats.mmseslope,agem_subject_stats.delta);
table1{7,1}='Number of subjects with baseline CDR=0';
table1{7,2}=sum(agem_subject_stats.baseline_cdr(agem_subject_stats.delta==0)==0);
table1{7,3}=sum(agem_subject_stats.baseline_cdr(agem_subject_stats.delta==1)==0);
table1{8,1}='Average length of time in study';
table1{8,2}=mean(agem_subject_stats.years_in_study(agem_subject_stats.delta==0));
table1{8,3}=mean(agem_subject_stats.years_in_study(agem_subject_stats.delta==1));
table1{8,4}=std(agem_subject_stats.years_in_study(agem_subject_stats.delta==0));
table1{8,5}=std(agem_subject_stats.years_in_study(agem_subject_stats.delta==1));
table1{8,6}=anova1(agem_subject_stats.years_in_study,agem_subject_stats.delta);
table1{9,1}='ApoE status';
table1{9,2}=sum(agem_subject_stats.apoe(agem_subject_stats.delta==0)==24) + sum(agem_subject_stats.apoe(agem_subject_stats.delta==0)==23) + sum(agem_subject_stats.apoe(agem_subject_stats.delta==0)==22);
table1{9,3}=sum(agem_subject_stats.apoe(agem_subject_stats.delta==1)==24) + sum(agem_subject_stats.apoe(agem_subject_stats.delta==1)==23) + sum(agem_subject_stats.apoe(agem_subject_stats.delta==1)==22);
table1{9,4}=sum(agem_subject_stats.apoe(agem_subject_stats.delta==0)==33); 
table1{9,5}=sum(agem_subject_stats.apoe(agem_subject_stats.delta==1)==33);
table1{9,6}=sum(agem_subject_stats.apoe(agem_subject_stats.delta==0)==34) + sum(agem_subject_stats.apoe(agem_subject_stats.delta==0)==44);
table1{9,7}=sum(agem_subject_stats.apoe(agem_subject_stats.delta==1)==34) + sum(agem_subject_stats.apoe(agem_subject_stats.delta==1)==44);

table1{10,1}='Race';
table1{10,2}=sum(strcmp(agem_subject_stats.Race(agem_subject_stats.delta==0),'Caucasian'));
table1{10,3}=sum(strcmp(agem_subject_stats.Race(agem_subject_stats.delta==1),'Caucasian'));
table1{10,4}=sum(strcmp(agem_subject_stats.Race(agem_subject_stats.delta==0),'African American'));
table1{10,5}=sum(strcmp(agem_subject_stats.Race(agem_subject_stats.delta==1),'African American'));
table1{11,1}='Ethnicity';
table1{11,2}=sum(strcmp(agem_subject_stats.Ethnicity(agem_subject_stats.delta==0),'Non-Hispanic'));
table1{11,3}=sum(strcmp(agem_subject_stats.Ethnicity(agem_subject_stats.delta==1),'Non-Hispanic'));
table1{11,4}=sum(strcmp(agem_subject_stats.Ethnicity(agem_subject_stats.delta==0),'Hispanic'));
table1{11,5}=sum(strcmp(agem_subject_stats.Ethnicity(agem_subject_stats.delta==1),'Hispanic'));

%[tbl,chi2,p]=crosstab(agem_subject_stats.Ethnicity,agem_subject_stats.delta)

%% apply linear mixed model to data
lme=fitlme(agem_repeated_amyloid, 'PET_fSUVR_rsf_TOT_CORTMEAN ~ 1+time_in_y*delta + baseline_cdr + baseline_age + sex + tracer + time_in_y+(1+time_in_y|Subject)');
agem_repeated_amyloid.fitted_PET_fSUVR_rsf_TOT_CORTMEAN=fitted(lme);


%% ploting the relationship of  PET_fSUVR_rsf_TOT_CORTMEAN for each subject - predicted values of amyloid deposition

figure;
for i=1:length(agem_unique_repeated_amyloid.Subject)
   id=strcmp(agem_repeated_amyloid.Subject,agem_unique_repeated_amyloid.Subject{i});
   if agem_repeated_amyloid.delta(id,1)==0 
   plot(agem_repeated_amyloid.age(id),agem_repeated_amyloid.fitted_PET_fSUVR_rsf_TOT_CORTMEAN(id),'LineWidth',0.75,'Color',[0 0.4 1]);hold on
   elseif agem_repeated_amyloid.delta(id,1)==1
   plot(agem_repeated_amyloid.age(id),agem_repeated_amyloid.fitted_PET_fSUVR_rsf_TOT_CORTMEAN(id),'LineWidth',1.75,'Color','r');hold on %% maintains plot so that later plots don't clear this data
   end 
        clear d
end
%% ploting the relationship of subcortical mean for each subject - predicted values of amyloid deposition

for i=1:length(agem_repeated_amyloid.Subject)
    subcort_values=table2array([agem_repeated_amyloid(i,17:20),agem_repeated_amyloid(i,135),agem_repeated_amyloid(i,150),agem_repeated_amyloid(i,152),agem_repeated_amyloid(i,168),agem_repeated_amyloid(i,170)]);
    agem_repeated_amyloid.PET_fSUVR_rsf_TOT_SUBCORTMEAN(i)=mean(subcort_values);
end

lme_subcort=fitlme(agem_repeated_amyloid, 'PET_fSUVR_rsf_TOT_SUBCORTMEAN ~ 1+time_in_y*delta + baseline_cdr + baseline_age + sex + tracer + time_in_y+(1+time_in_y|Subject)');
agem_repeated_amyloid.fitted_PET_fSUVR_rsf_TOT_SUBCORTMEAN=fitted(lme_subcort);

 figure;
 for i=1:length(agem_unique_repeated_amyloid.Subject)
    id=strcmp(agem_repeated_amyloid.Subject,agem_unique_repeated_amyloid.Subject{i});
    if agem_repeated_amyloid.delta(id,1)==0 
        plot(agem_repeated_amyloid.age(id),agem_repeated_amyloid.fitted_PET_fSUVR_rsf_TOT_SUBCORTMEAN(id),'LineWidth',0.75,'Color',[0 0.4 1]);hold on
    elseif agem_repeated_amyloid.delta(id,1)==1
        plot(agem_repeated_amyloid.age(id),agem_repeated_amyloid.fitted_PET_fSUVR_rsf_TOT_SUBCORTMEAN(id),'LineWidth',1.75,'Color','r');hold on
    end
    clear d
 end
 
figure;
plotResiduals(lme,'fitted');


[r1,r2,rEffects] = randomEffects(lme);
figure;
scatter(rEffects.Estimate(1:2:end),rEffects.Estimate(2:2:end))
title('Random Effects','FontSize',15)
xlabel('Intercept','FontSize',15)
ylabel('Slope','FontSize',15)

% slopes=rEffects.Estimate(2:2:end);


%% looping through all brain regions to see regional interactions 

agem_left_hemisphere=[agem_repeated_amyloid(:,2),agem_repeated_amyloid(:,5),agem_repeated_amyloid(:,7),agem_repeated_amyloid(:,277:283),agem_repeated_amyloid(:,43:45),agem_repeated_amyloid(:,47:48),agem_repeated_amyloid(:,50:52),agem_repeated_amyloid(:,54:59),agem_repeated_amyloid(:,61),agem_repeated_amyloid(:,60),agem_repeated_amyloid(:,62:75),agem_repeated_amyloid(:,49),agem_repeated_amyloid(:,76:77),agem_repeated_amyloid(:,53)];
agem_right_hemisphere=[agem_repeated_amyloid(:,2),agem_repeated_amyloid(:,5),agem_repeated_amyloid(:,7),agem_repeated_amyloid(:,277:283),agem_repeated_amyloid(:,91:93),agem_repeated_amyloid(:,95:96),agem_repeated_amyloid(:,98:100),agem_repeated_amyloid(:,102:107),agem_repeated_amyloid(:,109),agem_repeated_amyloid(:,108),agem_repeated_amyloid(:,110:123),agem_repeated_amyloid(:,97),agem_repeated_amyloid(:,124:125),agem_repeated_amyloid(:,101)];

agem_stats_table=table('Size',[68 7],'VariableTypes',{'string','double','double','double','double','double','double'});
agem_stats_table.Properties.VariableNames={'region','tStat_delta','Pval_delta','tStat_delta_interaction','Pval_delta_interaction','tStat_sex_interaction','Pval_sex_interaction'};

for i=11:44
    left_roi=[agem_left_hemisphere(:,1:10),agem_left_hemisphere(:,i)];
    left_roi.Properties.VariableNames={'Subject','tracer','time_in_y','delta','baseline_age','baseline_cdr','age','sex','education','race','amyloid_region'};
    lme=fitlme(left_roi, 'amyloid_region ~ 1+time_in_y*delta + 1+baseline_age*sex + baseline_cdr + baseline_age + sex + tracer + time_in_y+(1+time_in_y|Subject)');    
    agem_fitted_amyloid_regions(:,i-10)=fitted(lme);
    agem_stats_table.region(i-10,1)=agem_left_hemisphere.Properties.VariableNames(i);
    agem_stats_table.tStat_delta(i-10,1)=lme.Coefficients.tStat(4);
    agem_stats_table.Pval_delta(i-10,1)=lme.Coefficients.pValue(4);
    agem_stats_table.tStat_delta_interaction(i-10,1)=lme.Coefficients.tStat(8);
    agem_stats_table.Pval_delta_interaction(i-10,1)=lme.Coefficients.pValue(8);
    agem_stats_table.tStat_sex_interaction(i-10,1)=lme.Coefficients.tStat(9);
    agem_stats_table.Pval_sex_interaction(i-10,1)=lme.Coefficients.pValue(9);
end
for i=11:44
    right_roi=[agem_right_hemisphere(:,1:10),agem_right_hemisphere(:,i)];
    right_roi.Properties.VariableNames={'Subject','tracer','time_in_y','delta','baseline_age','baseline_cdr','age','sex','education','race','amyloid_region'};
    lme=fitlme(right_roi, 'amyloid_region ~ 1+time_in_y*delta + 1+baseline_age*sex + baseline_cdr + baseline_age + sex + tracer + time_in_y+(1+time_in_y|Subject)');    
    agem_fitted_amyloid_regions(:,i+24)=fitted(lme);
    agem_stats_table.region(i+24,1)=agem_right_hemisphere.Properties.VariableNames(i);
    agem_stats_table.tStat_delta(i+24,1)=lme.Coefficients.tStat(4);
    agem_stats_table.Pval_delta(i+24,1)=lme.Coefficients.pValue(4);
    agem_stats_table.tStat_delta_interaction(i+24,1)=lme.Coefficients.tStat(8);
    agem_stats_table.Pval_delta_interaction(i+24,1)=lme.Coefficients.pValue(8);
    agem_stats_table.tStat_sex_interaction(i+24,1)=lme.Coefficients.tStat(9);
    agem_stats_table.Pval_sex_interaction(i+24,1)=lme.Coefficients.pValue(9);
end

agem_subcortical_regions=[agem_repeated_amyloid(:,2),agem_repeated_amyloid(:,5),agem_repeated_amyloid(:,7),agem_repeated_amyloid(:,277:283),agem_repeated_amyloid(:,17:20),agem_repeated_amyloid(:,135),agem_repeated_amyloid(:,150),agem_repeated_amyloid(:,152),agem_repeated_amyloid(:,168),agem_repeated_amyloid(:,170)];
 for i=11:19
    sub_roi=[agem_subcortical_regions(:,1:10),agem_subcortical_regions(:,i)];
    sub_roi.Properties.VariableNames={'Subject','tracer','time_in_y','delta','baseline_age','baseline_cdr','age','sex','education','race','amyloid_region'};
    lme=fitlme(sub_roi, 'amyloid_region ~ 1+time_in_y*delta + 1+baseline_age*sex + baseline_cdr + baseline_age + sex + tracer + time_in_y+(1+time_in_y|Subject)');    
    agem_fitted_amyloid_regions(:,i+58)=fitted(lme);
    agem_stats_table.region(i+58,1)=agem_subcortical_regions.Properties.VariableNames(i);
    agem_stats_table.tStat_delta(i+58,1)=lme.Coefficients.tStat(4);
    agem_stats_table.Pval_delta(i+58,1)=lme.Coefficients.pValue(4);
    agem_stats_table.tStat_delta_interaction(i+58,1)=lme.Coefficients.tStat(8);
    agem_stats_table.Pval_delta_interaction(i+58,1)=lme.Coefficients.pValue(8);
    agem_stats_table.tStat_sex_interaction(i+58,1)=lme.Coefficients.tStat(9);1
    agem_stats_table.Pval_sex_interaction(i+58,1)=lme.Coefficients.pValue(9);
end
agem_stats_table.Pval_delta_fdr=mafdr(agem_stats_table.Pval_delta);
agem_stats_table.Pval_delta_inter_fdr=mafdr(agem_stats_table.Pval_delta_interaction);

%% only plot the significant regions
agem_stats_table.tStat_delta_inter_thresh=agem_stats_table.tStat_delta_interaction;
agem_stats_table.tStat_delta_inter_thresh(agem_stats_table.Pval_delta_inter_fdr>0.05)=0;

llabel_names=agem_stats_table.region(1:34);
rlabel_names=agem_stats_table.region(35:68);

%% find colour ranges for brain visualization map

% left side
clear full_cell
id=[1:35]'; 
for roi=1:length(llabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if agem_stats_table.tStat_delta_inter_thresh(roi)>=0
        blue_range=0:255;
        green_range=178.5:0.3:255;
    
    tstat_relative=1-agem_stats_table.tStat_delta_inter_thresh(roi)/4.235;
    color_relative=round(tstat_relative*256);
    
    R(roi,1)=255;
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=blue_range(color_relative);
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    elseif agem_stats_table.tStat_delta_inter_thresh(roi)<0;
        red_range=0:255;
        green_range=183.5:0.28:255;
    tstat_relative=1-agem_stats_table.tStat_delta_inter_thresh(roi)/-1.26;
    color_relative=round(tstat_relative*256);
    
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=255;
    A(roi,1)=0;
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end

% right side
clear full_cell
id=[1:35]'; 
for roi=1:length(rlabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if agem_stats_table.tStat_delta_inter_thresh(roi+34)>=0
        blue_range=0:255;
        green_range=178.5:0.3:255;
    
    tstat_relative=1-agem_stats_table.tStat_delta_inter_thresh(roi+34)/4.235;
    color_relative=round(tstat_relative*256);
    
    R(roi,1)=255;
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=blue_range(color_relative);
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=rlabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    elseif agem_stats_table.tStat_delta_inter_thresh(roi+34)<0;
        red_range=0:255;
        green_range=183.5:0.28:255;
    tstat_relative=1-agem_stats_table.tStat_delta_inter_thresh(roi+34)/-1.26;
    color_relative=round(tstat_relative*256);
    
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=255;
    A(roi,1)=0;
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=rlabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end
