%% amyloid slopes
agem_repeated_amyloid=sortrows(agem_repeated_amyloid,[2 6]);
agem_unique_repeated_amyloid=sortrows(agem_unique_repeated_amyloid,[2 6]);

% load('correlation_variables.mat');

unfitted_cortical_regions=[agem_repeated_amyloid(:,2),agem_repeated_amyloid(:,7),agem_repeated_amyloid(:,277:282),agem_repeated_amyloid(:,5),agem_repeated_amyloid(:,43:45),agem_repeated_amyloid(:,47:48),agem_repeated_amyloid(:,50:52),agem_repeated_amyloid(:,54:59),agem_repeated_amyloid(:,61),agem_repeated_amyloid(:,60),agem_repeated_amyloid(:,62:75),agem_repeated_amyloid(:,49),agem_repeated_amyloid(:,76:77),agem_repeated_amyloid(:,53),agem_repeated_amyloid(:,91:93),agem_repeated_amyloid(:,95:96),agem_repeated_amyloid(:,98:100),agem_repeated_amyloid(:,102:107),agem_repeated_amyloid(:,109),agem_repeated_amyloid(:,108),agem_repeated_amyloid(:,110:123),agem_repeated_amyloid(:,97),agem_repeated_amyloid(:,124:125),agem_repeated_amyloid(:,101),agem_repeated_amyloid(:,20)];
clear lme_table; lme_table=unfitted_cortical_regions(:,1:9);
fitted_cortical_regions=unfitted_cortical_regions(:,1:9);

for i=10:78
    lme_table.region=unfitted_cortical_regions{:,i};
    lme=fitlme(lme_table, 'region ~ 1+time_in_y*delta + 1+baseline_age*sex + baseline_cdr + baseline_age + sex + tracer + time_in_y+(1+time_in_y|Subject)');    
    fitted_cortical_regions{:,i}=fitted(lme);
    fitted_cortical_regions.Properties.VariableNames{i}=unfitted_cortical_regions.Properties.VariableNames{i};
end

amyloid_slopes=agem_unique_repeated_amyloid.Subject;
amyloid_slopes=cell2table(amyloid_slopes);
amyloid_slopes.Properties.VariableNames{1}='Subject';
clear C; clear ia; clear ic;
[C,ia,ic] = unique(fitted_cortical_regions.Subject,'first'); 
amyloid_parcellations_baseline = fitted_cortical_regions(ia,:);
clear C; clear ia; clear ic;
[C,ia,ic] = unique(fitted_cortical_regions.Subject,'last'); 
amyloid_parcellations_end = fitted_cortical_regions(ia,:); 

for i=10:78
    for j=1:length(agem_unique_repeated_amyloid.Subject)
        slope=(amyloid_parcellations_end{j,i}-amyloid_parcellations_baseline{j,i})/(amyloid_parcellations_end{j,2}-amyloid_parcellations_baseline{j,2});
        amyloid_slopes{j,i-8}=slope;
    end
    amyloid_slopes.Properties.VariableNames{i-8}=amyloid_parcellations_baseline.Properties.VariableNames{i};
end

amyloid_slopes.baseline_age=agem_unique_repeated_amyloid.baseline_age;
amyloid_slopes.delta=agem_unique_repeated_amyloid.delta;


thresh=abs(amyloid_slopes.PET_fSUVR_rsf_TOT_CAUD)<1;
figure; hold off;
    ix=amyloid_slopes.delta==0 & thresh==1; 
    l=0.5; k=1.5; orange=[1 0.5 0]; %https://medium.com/@SciencelyYours/matlab-colors-3ca3aa47040c
    y=amyloid_slopes.PET_fSUVR_rsf_TOT_CAUD(ix);
    x=amyloid_slopes.baseline_age(ix);
    [p,S]=polyfit(x,y,2); p_all(:,1)=p;  [Y,DELTA] = polyconf(p,x,S, 'alpha',0.1, 'predopt', 'curve');
    f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    fplot(f,'b', 'LineWidth',k);hold on;xlim([45, 95]); %scatter(age,MS(:,roi));  %xlim([20, 70]); ylim([-0.05 0.035]);     tmp=sortrows([x, Y+DELTA, Y-DELTA]); plot(tmp(:,1), tmp(:,2), 'b--', 'LineWidth',l);plot(tmp(:,1), tmp(:,3), 'b--', 'LineWidth',l);
    scatter(x,y, 10,'b', 'filled');
    x=[45:round(max(amyloid_slopes.baseline_age(ix)))]; x2=[x, fliplr(x)];[Y,DELTA] = polyconf(p,x,S, 'alpha',0.1, 'predopt', 'curve');
    inBetween = [Y-DELTA, fliplr(Y+DELTA)];
    fl1=fill(x2, inBetween, 'b'); alpha(fl1,.1);
    
     
    clear ix;
    ix=amyloid_slopes.delta==1 & thresh==1 
    l=0.5; k=1.5; orange=[1 0.5 0]; %https://medium.com/@SciencelyYours/matlab-colors-3ca3aa47040c
    y=amyloid_slopes.PET_fSUVR_rsf_TOT_CAUD(ix);
    x=amyloid_slopes.baseline_age(ix);
    [p,S]=polyfit(x,y,2); p_all(:,1)=p;  [Y,DELTA] = polyconf(p,x,S, 'alpha',0.1, 'predopt', 'curve');
    f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    fplot(f,'r', 'LineWidth',k);hold on;xlim([45, 95]); %scatter(age,MS(:,roi));  %xlim([20, 70]); ylim([-0.05 0.035]);
    tmp=sortrows([x, Y+DELTA, Y-DELTA]); plot(tmp(:,1), tmp(:,2), 'r--', 'LineWidth',l);plot(tmp(:,1), tmp(:,3), 'r--', 'LineWidth',l); hold on;
    scatter(x,y, 10,'r', 'filled');hold on;
    x=[45:round(max(amyloid_slopes.baseline_age(ix)))]; x2=[x, fliplr(x)];[Y,DELTA] = polyconf(p,x,S, 'alpha',0.1, 'predopt', 'curve');
    inBetween = [Y-DELTA, fliplr(Y+DELTA)];
    fl1=fill(x2, inBetween, 'r'); alpha(fl1,.1);
    
   

%% fmri slopes
fmri_metrics_run1=sortrows(fmri_metrics_run1,[1 3]);
unique_fmri_subs=sortrows(unique_fmri_subs,[1 3]);


clear lme_table; lme_table=[fmri_metrics_run1(:,1),fmri_metrics_run1(:,3),fmri_metrics_run1(:,13:14),fmri_metrics_run1(:,437:444)];
fitted_connections=lme_table(:,1:2);
for i=17:226
    lme_table.connect=fmri_metrics_run1{:,i};
    lme=fitlme(lme_table, 'connect ~ 1+time_in_y*delta + prop_signal + mean_fd + baseline_cdr + baseline_age + sex + time_in_y + (1+time_in_y|Subs)');
    fitted_connections{:,i-14}=fitted(lme);
    fitted_connections.Properties.VariableNames{i-14}=fmri_metrics_run1.Properties.VariableNames{i};
end

fmri_slopes=unique_fmri_subs.Subs;
fmri_slopes=sortrows(fmri_slopes,[1]);
fmri_slopes=cell2table(fmri_slopes);
fmri_slopes.Properties.VariableNames{1}='Subject';
clear C; clear ia; clear ic;
[C,ia,ic] = unique(fitted_connections.Subs,'first'); 
fmri_parcellations_baseline = fitted_connections(ia,:); 
clear C; clear ia; clear ic;
[C,ia,ic] = unique(fitted_connections.Subs,'last'); 
fmri_parcellations_end = fitted_connections(ia,:); 

for i=3:212
    for j=1:length(fmri_parcellations_end.Subs)
        slope=(fmri_parcellations_end{j,i}-fmri_parcellations_baseline{j,i})/((fmri_parcellations_end{j,2}-fmri_parcellations_baseline{j,2})/365);
        fmri_slopes{j,i-1}=slope;
    end
    fmri_slopes.Properties.VariableNames{i-1}=fmri_parcellations_baseline.Properties.VariableNames{i};
end

fmri_slopes.baseline_age=unique_fmri_subs.baseline_age;
fmri_slopes.delta=unique_fmri_subs.delta;

thresh=abs(fmri_slopes.F_IC5IC11)<3;
figure(7); hold off;
    clear ix;
    ix=fmri_slopes.delta==0 & thresh==1; 
    l=0.5; k=1.5; orange=[1 0.5 0]; %https://medium.com/@SciencelyYours/matlab-colors-3ca3aa47040c
    y=fmri_slopes.F_IC5IC11(ix);
    x=fmri_slopes.baseline_age(ix);
    [p,S]=polyfit(x,y,2); p_all(:,1)=p;  [Y,DELTA] = polyconf(p,x,S, 'alpha',0.1, 'predopt', 'curve');
    f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    figure(7);fplot(f,'b', 'LineWidth',k);hold on;xlim([45, 95]); %scatter(age,MS(:,roi));  %xlim([20, 70]); ylim([-0.05 0.035]);     tmp=sortrows([x, Y+DELTA, Y-DELTA]); plot(tmp(:,1), tmp(:,2), 'b--', 'LineWidth',l);plot(tmp(:,1), tmp(:,3), 'b--', 'LineWidth',l);
    scatter(x,y, 10,'b', 'filled');
    x=[45:round(max(fmri_slopes.baseline_age(ix)))]; x2=[x, fliplr(x)];[Y,DELTA] = polyconf(p,x,S, 'alpha',0.1, 'predopt', 'curve');
    inBetween = [Y-DELTA, fliplr(Y+DELTA)];
    fl1=fill(x2, inBetween, 'b'); alpha(fl1,.1);

    clear ix;
    ix=fmri_slopes.delta==1 & thresh==1; 
    l=0.5; k=1.5; orange=[1 0.5 0]; %https://medium.com/@SciencelyYours/matlab-colors-3ca3aa47040c
    y=fmri_slopes.F_IC5IC11(ix);
    x=fmri_slopes.baseline_age(ix);
    [p,S]=polyfit(x,y,2); p_all(:,1)=p;  [Y,DELTA] = polyconf(p,x,S, 'alpha',0.1, 'predopt', 'curve');
    f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    figure(7); fplot(f,'r', 'LineWidth',k);hold on;xlim([45, 95]); %scatter(age,MS(:,roi));  %xlim([20, 70]); ylim([-0.05 0.035]);
    tmp=sortrows([x, Y+DELTA, Y-DELTA]); plot(tmp(:,1), tmp(:,2), 'r--', 'LineWidth',l);plot(tmp(:,1), tmp(:,3), 'r--', 'LineWidth',l); hold on;
    scatter(x,y, 10,'r', 'filled');hold on;
    x=[45:round(max(fmri_slopes.baseline_age(ix)))]; x2=[x, fliplr(x)];[Y,DELTA] = polyconf(p,x,S, 'alpha',0.1, 'predopt', 'curve');
    inBetween = [Y-DELTA, fliplr(Y+DELTA)];
    fl1=fill(x2, inBetween, 'r'); alpha(fl1,.1);
     
    
    



%% correlation

matched_amyloid_slopes=amyloid_slopes(:,:);
matched_amyloid_slopes(~ismember(amyloid_slopes.Subject,fmri_slopes.Subject),:)=[];
matched_slopes=[matched_amyloid_slopes(:,1:72),fmri_slopes(:,2:211)];

corr([matched_slopes.PET_fSUVR_rsf_L_CTX_FUSIFORM,matched_slopes.PET_fSUVR_rsf_R_CTX_FUSIFORM,matched_slopes.PET_fSUVR_rsf_L_CTX_INFRPRTL,matched_slopes.PET_fSUVR_rsf_R_CTX_INFPRTL,matched_slopes.PET_fSUVR_rsf_L_CTX_INFRTMP,matched_slopes.PET_fSUVR_rsf_R_CTX_INFTMP,matched_slopes.PET_fSUVR_rsf_L_CTX_LATOCC,matched_slopes.PET_fSUVR_rsf_R_CTX_LATOCC,matched_slopes.PET_fSUVR_rsf_L_CTX_MIDTMP,matched_slopes.PET_fSUVR_rsf_R_CTX_MIDTMP,matched_slopes.PET_fSUVR_rsf_L_CTX_PARACNTRL,matched_slopes.PET_fSUVR_rsf_R_CTX_PARACNTRL,matched_slopes.PET_fSUVR_rsf_R_CTX_PARSTRNGLRS,matched_slopes.PET_fSUVR_rsf_L_CTX_POSTCNTRL,matched_slopes.PET_fSUVR_rsf_R_CTX_POSTCNTRL,matched_slopes.PET_fSUVR_rsf_R_CTX_PRECNTRL,matched_slopes.PET_fSUVR_rsf_L_CTX_SUPERTMP,matched_slopes.PET_fSUVR_rsf_R_CTX_SUPERTMP]);
[coeff,scores,latent,~,explained]= pca([matched_slopes.PET_fSUVR_rsf_L_CTX_FUSIFORM,matched_slopes.PET_fSUVR_rsf_R_CTX_FUSIFORM,matched_slopes.PET_fSUVR_rsf_L_CTX_INFRPRTL,matched_slopes.PET_fSUVR_rsf_R_CTX_INFPRTL,matched_slopes.PET_fSUVR_rsf_L_CTX_INFRTMP,matched_slopes.PET_fSUVR_rsf_R_CTX_INFTMP,matched_slopes.PET_fSUVR_rsf_L_CTX_LATOCC,matched_slopes.PET_fSUVR_rsf_R_CTX_LATOCC,matched_slopes.PET_fSUVR_rsf_L_CTX_MIDTMP,matched_slopes.PET_fSUVR_rsf_R_CTX_MIDTMP,matched_slopes.PET_fSUVR_rsf_L_CTX_PARACNTRL,matched_slopes.PET_fSUVR_rsf_R_CTX_PARACNTRL,matched_slopes.PET_fSUVR_rsf_R_CTX_PARSTRNGLRS,matched_slopes.PET_fSUVR_rsf_L_CTX_POSTCNTRL,matched_slopes.PET_fSUVR_rsf_R_CTX_POSTCNTRL,matched_slopes.PET_fSUVR_rsf_R_CTX_PRECNTRL,matched_slopes.PET_fSUVR_rsf_L_CTX_SUPERTMP,matched_slopes.PET_fSUVR_rsf_R_CTX_SUPERTMP]);
matched_slopes.scores1=scores(:,1);    

figure();
clear ix; ix = matched_slopes.delta==0 & abs(matched_slopes.F_IC5IC11)<2 & abs(matched_slopes.scores1)<10;
scatter(matched_slopes.F_IC6IC21(ix),matched_slopes.scores1(ix),10,'filled','b');hold on;
clear ix; ix = matched_slopes.delta==1 & abs(matched_slopes.F_IC6IC21)<2 & abs(matched_slopes.scores1)<10;
scatter(matched_slopes.F_IC6IC21(ix),matched_slopes.scores1(ix),10,'filled','r');hold on;
clear ix; ix=abs(matched_slopes.F_IC2IC6)<2 & abs(matched_slopes.scores1)<10;
coef=polyfit(matched_slopes.F_IC6IC21(ix),matched_slopes.scores1(ix),1);
refline(coef(1),coef(2));


clear ix;ix=abs(matched_slopes.F_IC14IC21)<2 & abs(matched_slopes.scores1)<1;
[r,p]=corrcoef (matched_slopes.F_IC14IC21(ix),matched_slopes.scores1(ix)) %% repeat for all sign connectivities
sum(ix) % -2 = degrees of freedom

corrtable=readtable('OASIS3_papertables.xlsx','Sheet','Correlation_tables');
corrtable.fdr_PValue=mafdr(corrtable.PValue,'BHFDR','true');
