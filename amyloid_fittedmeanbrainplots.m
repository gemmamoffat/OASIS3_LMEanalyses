agem_fitted_amyloid_regions=array2table(agem_fitted_amyloid_regions);
agem_fitted_amyloid_regions(:,78)=agem_left_hemisphere(:,1);
agem_fitted_amyloid_regions.Properties.VariableNames{78}='Subject';
agem_fitted_amyloid_regions(:,79)=agem_left_hemisphere(:,4);
agem_fitted_amyloid_regions.Properties.VariableNames{79}='delta';
% fitted_amyloid_regions(:,71)=repeated_amyloid(:,291);
% fitted_amyloid_regions.Properties.VariableNames{71}='delta_tertiary';
% add column for secondary delta

clear C; clear ia; clear ic;
[C,ia,ic] = unique(agem_fitted_amyloid_regions(:,78),'first'); 
baseline_fitted_amyloid_regions = agem_fitted_amyloid_regions(ia,:);
baseline0_regions=baseline_fitted_amyloid_regions(baseline_fitted_amyloid_regions.delta==0,:);
baseline1_regions=baseline_fitted_amyloid_regions(baseline_fitted_amyloid_regions.delta==1,:);
clear C; clear ia; clear ic;
[C,ia,ic] = unique(agem_fitted_amyloid_regions(:,78),'last'); 
end_fitted_amyloid_regions = agem_fitted_amyloid_regions(ia,:);
end0_regions=end_fitted_amyloid_regions(end_fitted_amyloid_regions.delta==0,:);
end1_regions=end_fitted_amyloid_regions(end_fitted_amyloid_regions.delta==1,:);

baseline0_stats=agem_stats_table(:,{'region'});
baseline1_stats=agem_stats_table(:,{'region'});
end0_stats=agem_stats_table(:,{'region'});
end1_stats=agem_stats_table(:,{'region'});

for i=1:77
    baseline0_stats.mean_fitted_value(i)=mean(baseline0_regions{:,i});
    baseline1_stats.mean_fitted_value(i)=mean(baseline1_regions{:,i});
    end0_stats.mean_fitted_value(i)=mean(end0_regions{:,i});
    end1_stats.mean_fitted_value(i)=mean(end1_regions{:,i});
end

%% baseline non-progressor values;, left side
clear full_cell
id=[1:35]'; 
for roi=1:length(llabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if baseline0_stats.mean_fitted_value(roi)>=0
        blue_range=102:0.75:255;
        green_range=51:0.87:229;
        red_range=0:204;
    mean_relative=1-(baseline0_stats.mean_fitted_value(roi)-0.9)/(2.3-0.9);
    color_relative=round(mean_relative*205);
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=round(blue_range(color_relative));
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end
%% baseline non-progressor values; right side
clear full_cell
id=[1:35]'; 
for roi=1:length(rlabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if baseline0_stats.mean_fitted_value(roi+34)>=0
        blue_range=102:0.75:255;
        green_range=51:0.87:229;
        red_range=0:204;
    mean_relative=1-(baseline0_stats.mean_fitted_value(roi+34)-0.9)/(2.3-0.9);
    color_relative=round(mean_relative*205);
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=round(blue_range(color_relative));
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end

%% baseline progressor values, left side
clear full_cell
id=[1:34]'; 
for roi=1:length(llabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if baseline1_stats.mean_fitted_value(roi)>=0
        blue_range=102:0.75:255;
        green_range=51:0.87:229;
        red_range=0:204;
    mean_relative=1-(baseline1_stats.mean_fitted_value(roi)-0.9)/(2.3-0.9);
    color_relative=round(mean_relative*205);
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=round(blue_range(color_relative));
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end

%% baseline progressors; right side
clear full_cell
id=[1:35]'; 
for roi=1:length(rlabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if baseline1_stats.mean_fitted_value(roi+34)>=0
        blue_range=102:0.75:255;
        green_range=51:0.87:229;
        red_range=0:204;
    mean_relative=1-(baseline1_stats.mean_fitted_value(roi+34)-0.9)/(2.3-0.9);
    color_relative=round(mean_relative*205);
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=round(blue_range(color_relative));
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=rlabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end


%% end non-progressor values; left side
clear full_cell
id=[1:34]'; 
for roi=1:length(llabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if end0_stats.mean_fitted_value(roi)>=0
        blue_range=102:0.75:255;
        green_range=51:0.87:229;
        red_range=0:204;
    mean_relative=1-(end0_stats.mean_fitted_value(roi)-0.9)/(2.3-0.9);
    color_relative=round(mean_relative*205);
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=round(blue_range(color_relative));
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end

%% end non-progressor values; right side
clear full_cell
id=[1:35]'; 
for roi=1:length(rlabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if end0_stats.mean_fitted_value(roi+34)>=0
        blue_range=102:0.75:255;
        green_range=51:0.87:229;
        red_range=0:204;
    mean_relative=1-(end0_stats.mean_fitted_value(roi+34)-0.9)/(2.3-0.9);
    color_relative=round(mean_relative*205);
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=round(blue_range(color_relative));
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=rlabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end

%% end progressor values; left side
clear full_cell
id=[1:34]'; 
for roi=1:length(llabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if end1_stats.mean_fitted_value(roi)>=0
        blue_range=102:0.75:255;
        green_range=51:0.87:229;
        red_range=0:204;
    mean_relative=1-(end1_stats.mean_fitted_value(roi)-0.9)/(2.3-0.9);
    color_relative=round(mean_relative*205);
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=round(blue_range(color_relative));
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=llabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end

%% end progressor values; right side
clear full_cell
id=[1:35]'; 
for roi=1:length(rlabel_names)
    %e.g. t=2; maxt=4 => then color should be half way to the max of the range:
    if end1_stats.mean_fitted_value(roi+34)>=0
        blue_range=102:0.75:255;
        green_range=51:0.87:229;
        red_range=0:204;
    mean_relative=1-(end1_stats.mean_fitted_value(roi+34)-0.9)/(2.3-0.9);
    color_relative=round(mean_relative*205);
    
    R(roi,1)=red_range(color_relative);
    G(roi,1)=round(green_range(color_relative));
    B(roi,1)=round(blue_range(color_relative));
    A(roi,1)=0;
    
    full_cell{roi, 1}=id(roi); full_cell{roi, 2}=rlabel_names{roi};full_cell{roi, 3}=R(roi);    full_cell{roi, 4}=G(roi);    full_cell{roi, 5}=B(roi);     full_cell{roi, 6}=A(roi);
    end
end