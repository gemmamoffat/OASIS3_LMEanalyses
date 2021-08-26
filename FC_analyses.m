
fmri_metrics_run1=readtable('fmri_preproc_data_subsettedrun1.xlsx');
%load('/Users/gemmamoffat/Desktop/CAMH/reports/amyloid_variables_variables.mat')

for i=1:length(fmri_metrics_run1.Subs)
    if ismember(fmri_metrics_run1.Subs{i},allsubs)
        id=strcmp(fmri_metrics_run1.Subs{i},allsubs);
        fmri_metrics_run1.delta(i,1)=delta(id==1);
        fmri_metrics_run1.baseline_age(i,1)=age(id==1);
        fmri_metrics_run1.baseline_cdr(i,1)=baseline_cdr(id==1);
        fmri_metrics_run1.age(i,1)=fmri_metrics_run1.Ses(i)/365 +age(id==1);
        fmri_metrics_run1.time_in_y(i,1)=fmri_metrics_run1.Ses(i)/365;
        fmri_metrics_run1.sex(i,1)=sex(id==1);
        fmri_metrics_run1.education(i,1)=education(id==1);
        fmri_metrics_run1.race(i,1)=race(id==1);
    end
end


fmri_metrics_run1(~ismember(fmri_metrics_run1.Subs,agem_unique_repeated_amyloid.Subject),:)=[];

fmri_metrics_run1=grouptransform(fmri_metrics_run1,'Subs',@numel,'ReplaceValues', false); % counts how many measures per subject, adds column of vaues
fmri_metrics_run1=fmri_metrics_run1(1:length(fmri_metrics_run1.Subs),1:445); % only keeps valid columns
fmri_metrics_run1.Properties.VariableNames{445} = 'num_sessions';
fmri_metrics_run1(fmri_metrics_run1.num_sessions==1 | fmri_metrics_run1.mean_fd>0.1,:)=[];

[C,ia,ic] = unique(fmri_metrics_run1.Subs); 
unique_fmri_subs=fmri_metrics_run1(ia,:);

fc_stats_table=table('Size',[210 4],'VariableTypes',{'double','double','double','double'});
fc_stats_table.Properties.VariableNames={'tStat_delta','Pval_delta','tStat_delta_interaction','Pval_delta_interaction'};

for i=17:226
    lme_table=[fmri_metrics_run1(:,1),fmri_metrics_run1(:,13:14),fmri_metrics_run1(:,437:445)];
    lme_table.connection=fmri_metrics_run1{:,i};
    lme=fitlme(lme_table, 'connection ~ 1+time_in_y*delta + prop_signal + mean_fd + baseline_cdr + baseline_age + sex + time_in_y + (1+time_in_y|Subs)');
    fc_stats_table.tStat_delta(i-16,1)=lme.Coefficients.tStat(4);
    fc_stats_table.Pval_delta(i-16,1)=lme.Coefficients.pValue(4);
    fc_stats_table.tStat_delta_interaction(i-16,1)=lme.Coefficients.tStat(9);
    fc_stats_table.Pval_delta_interaction(i-16,1)=lme.Coefficients.pValue(9);
end

lme_info = [fmri_metrics_run1(:,1),fmri_metrics_run1(:,13:14),fmri_metrics_run1(:,437:445),fmri_metrics_run1(:,17:226)];
% instead of using mafdr; moved data to R Studio and completed permutation
% testing on it
% reimport data back to MATLAB to visualize
permutted_stats=readtable('OASIS3_papertables.xlsx','Sheet','FC_interaction');
for i=1:length(fc_stats_table.Pval_delta_interaction)
    if fc_stats_table.Pval_delta_interaction(i) < 0.02
        disp(i);
    end
end
tstats=fc_stats_table.tStat_delta;
tstats(fc_stats_table.Pval_delta>=0.025) = 0;
tstats(12,1)=0;
tstats(20,1)=0;
tstats(94,1)=0;
tstats(99,1)=0;
tstats(154,1)=0;
tstats(155,1)=0;
tstats(159,1)=0;


%% 
square_net_from_file = zeros(21);
square_net_from_file(triu(ones(21),1)>0) = tstats;
ica2yeo7=readtable('ica2yeo7.csv');
figure; imagesc(square_net_from_file); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', ica2yeo7.Yeo7N, 'XTickLabelRotation',90);
set(gca, 'YTick', 1:21, 'YTickLabel', ica2yeo7.Yeo7N);


%% disclding some regions
% exclude_cr=[0	1	0	0	0	0	1	0	1	0	1	0	1	0	1	1	0	1	0	1	1	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	0	1	0	1	1	1	0	1	1	0	1	0	1	1	1	0	1	1	1	0	1	0	1	1	1	0	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
% exclude_cr=logical(exclude_cr');
% fc_stats_table_relevant=fc_stats_table(exclude_cr,:);
% 
% fc_stats_table_relevant.Pval_delta_fdr=mafdr(fc_stats_table_relevant.Pval_delta);
% fc_stats_table_relevant.Pval_delta_inter_fdr=mafdr(fc_stats_table_relevant.Pval_delta_interaction);
% for i=1:length(fc_stats_table_relevant.Pval_delta_interaction)
%     if fc_stats_table_relevant.Pval_delta_inter_fdr(i) < 0.1
%         disp(i);
%     end
% end
% tstats=fc_stats_table_relevant.tStat_delta_interaction;
% tstats(fc_stats_table_relevant.Pval_delta_inter_fdr>=0.1) = 0;
%% 
% square_net_from_file = zeros(15);
% square_net_from_file(triu(ones(15),1)>0) = tstats;
% ica2yeo7=readtable('ica2yeo7.csv');
% figure(3); imagesc(square_net_from_file); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', ica2yeo7.Yeo7N, 'XTickLabelRotation',90);
% set(gca, 'YTick', 1:21, 'YTickLabel', ica2yeo7.Yeo7N);
%https://www.fmrib.ox.ac.uk/ukbiobank/group_means/netmat_info.html

%% lme loop - partial 

% pc_stats_table=table('Size',[210 4],'VariableTypes',{'double','double','double','double'});
% pc_stats_table.Properties.VariableNames={'tStat_delta','Pval_delta','tStat_delta_interaction','Pval_delta_interaction'};
% 
% for i=227:436
%     lme_table=[fmri_metrics_run1(:,1),fmri_metrics_run1(:,13:14),fmri_metrics_run1(:,437:445)];
%     lme_table.connection=fmri_metrics_run1{:,i};
%     lme=fitlme(lme_table, 'connection ~ 1+time_in_y*delta + prop_signal + mean_fd + baseline_cdr + baseline_age + sex + time_in_y+(1+time_in_y|Subs)');
%     pc_stats_table.tStat_delta(i-226,1)=lme.Coefficients.tStat(4);
%     pc_stats_table.Pval_delta(i-226,1)=lme.Coefficients.pValue(4);
%     pc_stats_table.tStat_delta_interaction(i-226,1)=lme.Coefficients.tStat(9);
%     pc_stats_table.Pval_delta_interaction(i-226,1)=lme.Coefficients.pValue(9);
% end
% 
% exclude_cr=[0	1	0	0	0	0	1	0	1	0	1	0	1	0	1	1	0	1	0	1	1	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	0	1	0	1	1	1	0	1	1	0	1	0	1	1	1	0	1	1	1	0	1	0	1	1	1	0	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
% exclude_cr=logical(exclude_cr');
% pc_stats_table_relevant=pc_stats_table(exclude_cr,:);
% 
% pc_stats_table_relevant.Pval_delta_fdr=mafdr(pc_stats_table_relevant.Pval_delta);
% pc_stats_table_relevant.Pval_delta_inter_fdr=mafdr(pc_stats_table_relevant.Pval_delta_interaction);
% 
% for i=1:length(pc_stats_table_relevant.Pval_delta_interaction)
%     if pc_stats_table_relevant.Pval_delta_inter_fdr(i) < 0.1
%         disp(i);
%     end
% end
% 
% tstats=pc_stats_table_relevant.tStat_delta_interaction;
% tstats(pc_stats_table_relevant.Pval_delta_inter_fdr >=0.1) = 0;
% 
% square_net_from_file = zeros(15);
% square_net_from_file(triu(ones(15),1)>0) = tstats;
% figure(4); imagesc(square_net_from_file); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', ica2yeo7.Yeo7N, 'XTickLabelRotation',90);
% set(gca, 'YTick', 1:21, 'YTickLabel', ica2yeo7.Yeo7N);

%% spagetti plot
%  
    lme_table=[fmri_metrics_run1(:,1),fmri_metrics_run1(:,13:14),fmri_metrics_run1(:,437:445)];
    lme_table.connection=fmri_metrics_run1{:,66};
    lme=fitlme(lme_table, 'connection ~ 1+time_in_y*delta + prop_signal + mean_fd + baseline_cdr + baseline_age + sex + time_in_y+(1+time_in_y|Subs)');
    fmri_metrics_run1.fitted_connection50=fitted(lme);
    fmri_subs=unique(fmri_metrics_run1.Subs);
    
    figure;
    for i=1:length(unique_fmri_subs.Subs)
       id=strcmp(fmri_metrics_run1.Subs,fmri_subs{i});
       if fmri_metrics_run1.delta(id,1)==0 
       plot(fmri_metrics_run1.age(id),fmri_metrics_run1.fitted_connection50(id),'LineWidth',0.75,'Color',[0 0.4 1]);hold on
       elseif fmri_metrics_run1.delta(id,1)==1
       plot(fmri_metrics_run1.age(id),fmri_metrics_run1.fitted_connection50(id),'LineWidth',1.75,'Color','r');hold on %% maintains plot so that later plots don't clear this data
       end 
    end

 
% %% secondary analyses
% fc_stats_table2=table('Size',[210 4],'VariableTypes',{'double','double','double','double'});
% fc_stats_table2.Properties.VariableNames={'tStat_delta','Pval_delta','tStat_delta_interaction','Pval_delta_interaction'};
% 
% for i=17:226
%     lme_table=[fmri_metrics_run1(:,1),fmri_metrics_run1(:,13:14),fmri_metrics_run1(:,437:446)];
%     lme_table.connection=fmri_metrics_run1{:,i};
%     lme=fitlme(lme_table, 'connection ~ 1+time_in_y*delta_secondary + prop_signal + mean_fd + baseline_cdr + baseline_age + sex + time_in_y + (1+time_in_y|Subs)');
%     fc_stats_table2.tStat_delta(i-16,1)=lme.Coefficients.tStat(8);
%     fc_stats_table2.Pval_delta(i-16,1)=lme.Coefficients.pValue(8);
%     fc_stats_table2.tStat_delta_interaction(i-16,1)=lme.Coefficients.tStat(9);
%     fc_stats_table2.Pval_delta_interaction(i-16,1)=lme.Coefficients.pValue(9);
% end
% fc_stats_table2.Pval_delta_fdr=mafdr(fc_stats_table2.Pval_delta);
% fc_stats_table2.Pval_delta_inter_fdr=mafdr(fc_stats_table2.Pval_delta_interaction);
% 
% 
% for i=1:length(fc_stats_table2.Pval_delta_interaction)
%     if fc_stats_table2.Pval_delta(i) < 0.01
%         disp(i);
%     end
% end
% tstats=fc_stats_table2.tStat_delta;
% tstats(fc_stats_table2.Pval_delta>=0.01) = 0;
% %% 
% square_net_from_file = zeros(21);
% square_net_from_file(triu(ones(21),1)>0) = tstats;
% ica2yeo7=readtable('ica2yeo7.csv');
% figure(3); imagesc(square_net_from_file); colorbar; caxis([-5 5]); set(gca, 'XTick', 1:21, 'XTickLabel', ica2yeo7.Yeo7N, 'XTickLabelRotation',90);
% set(gca, 'YTick', 1:21, 'YTickLabel', ica2yeo7.Yeo7N);

