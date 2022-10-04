clearvars -except c awakeind_match locind_match
% clearvars -except c
animal = 'p9';
load(strcat(animal,'_V1_spikes_filters.mat'))
load(strcat(animal,'_LGN_spikes_filters.mat'))
load(strcat(animal,'_LM_spikes_filters.mat'))
load(strcat(animal,'_awake_info.mat'))
load(strcat(animal,'_loc_info.mat'))

ext = 'shuff';
win = 1;
nsub = 5;
un_oris = [0 30 60 90];

ee = stim_periods;

% locvel = [];
% for i = 1:size(locind,2)
%     [val,idx1]=min(abs(veltime - ee(locind(i),1)));
%     [val,idx2]=min(abs(veltime - ee(locind(i),2)));
%     locvel(i) = mean(v2(idx1:idx2));
% end
% f = find(locvel >= 2);
% locvel = locvel(f);
% locind = locind(f);

loc_phase_2 = ss2;
fano_cut = 4;
V_mat = V_mat((V_fanos <= fano_cut),:,:);
fano_cut = 4;
LGN_mat = LGN_mat((LGN_fanos <= fano_cut),:,:);
fano_cut = 5;
LM_mat = LM_mat((LM_fanos <= fano_cut),:,:);

oritemp = ori(locind);%assumes locomotion trials < stationary trials
s = [];
for i = 1:length(un_oris)
    s(i) = sum(oritemp == un_oris(i));
end
s
n_ori = min(s);

n = input('Enter a number, 1 for locomotion, 0 for stationary: '); 
switch n
    case 1
        ee2 = ee(locind,:);
        ori2 = ori(locind);
        V_mat = V_mat(:,locind,:);
        LGN_mat = LGN_mat(:,locind,:);
        LM_mat = LM_mat(:,locind,:);
        if exist('locind_match','var') == 0
            locind_match = [];
            for i = 1:length(un_oris)
                f = find(ori2 == un_oris(i));
                locind_match = [locind_match f(randperm(length(f),n_ori))];
            end
            locind_match = sort(locind_match);
        end
        ee2 = ee2(locind_match,:);
        ori2 = ori2(locind_match);
        V_mat = V_mat(:,locind_match,:);
        LGN_mat = LGN_mat(:,locind_match,:);
        LM_mat = LM_mat(:,locind_match,:);
    case 0
        ee2 = ee(awakeind,:);
        ori2 = ori(awakeind);
        V_mat = V_mat(:,awakeind,:);
        LGN_mat = LGN_mat(:,awakeind,:);
        LM_mat = LM_mat(:,awakeind,:);     
        if exist('awakeind_match','var') == 0
            awakeind_match = [];
            for i = 1:length(un_oris)
                f = find(ori2 == un_oris(i));
                awakeind_match = [awakeind_match f(randperm(length(f),n_ori))];
            end
            awakeind_match = sort(awakeind_match);
        end
        ee2 = ee2(awakeind_match,:);
        ori2 = ori2(awakeind_match);
        V_mat = V_mat(:,awakeind_match,:);
        LGN_mat = LGN_mat(:,awakeind_match,:);
        LM_mat = LM_mat(:,awakeind_match,:);
end

%%
winwidth = 150;
winstep = 10;
V_mat2 = [];
j = 1;
for i = 1:winstep:(size(V_mat,3)-winwidth)
    V_mat2(:,:,j) = sum(V_mat(:,:,i:i+(winwidth-1)),3);
    LGN_mat2(:,:,j) = sum(LGN_mat(:,:,i:i+(winwidth-1)),3);
    LM_mat2(:,:,j) = sum(LM_mat(:,:,i:i+(winwidth-1)),3);
    j = j+1;
end

% remove nonvarying neurons and trials
% for k = 1:size(V_mat2,3)
%     for i = 1:size(V_mat2,1)
%         if sum(V_mat2(i,:,k)) == 0
%             V_mat2(i,:,k) = NaN;
%         end
%     end        
%     for j = 1:size(V_mat2,2)
%         if sum(V_mat2(:,j,k)) == 0
%             V_mat2(:,j,k) = NaN;
%         end
%     end
% end

f_ori = [];
for j = 1:length(un_oris)
    f_ori = [f_ori find(ori2 == [un_oris(j)])];
end
f_ori = sort(f_ori);
ori3 = ori2(f_ori);

V_mat3 = V_mat2(:,f_ori,:);
LGN_mat3 = LGN_mat2(:,f_ori,:);
LM_mat3 = LM_mat2(:,f_ori,:);

V_matonez = zeros(size(V_mat3,1),size(V_mat3,2),size(V_mat3,3));
V_matshuff = zeros(size(V_mat3,1),size(V_mat3,2),size(V_mat3,3));
V_matallz = zeros(size(V_mat3,1),size(V_mat3,2),size(V_mat3,3));
LGN_matonez = zeros(size(LGN_mat3,1),size(LGN_mat3,2),size(LGN_mat3,3));
LGN_matshuff = zeros(size(LGN_mat3,1),size(LGN_mat3,2),size(LGN_mat3,3));
LGN_matallz = zeros(size(LGN_mat3,1),size(LGN_mat3,2),size(LGN_mat3,3));
LM_matonez = zeros(size(LM_mat3,1),size(LM_mat3,2),size(LM_mat3,3));
LM_matshuff = zeros(size(LM_mat3,1),size(LM_mat3,2),size(LM_mat3,3));
LM_matallz = zeros(size(LM_mat3,1),size(LM_mat3,2),size(LM_mat3,3));

for i = 1:size(V_mat3,1)
    for k = 1:size(V_mat3,3)
        V_matonez(i,:,k) = zscore(V_mat3(i,:,k));  
        for j = 1:length(un_oris)
            f_ori = find(ori3 == [un_oris(j)]);
            V_matallz(i,f_ori,k) = zscore(V_mat3(i,f_ori,k));  
            V_matshuff(i,f_ori,k) = V_matallz(i,f_ori(randperm(length(f_ori))),k);
        end
    end
end
for i = 1:size(LGN_mat3,1)
    for k = 1:size(LGN_mat3,3)
        LGN_matonez(i,:,k) = zscore(LGN_mat3(i,:,k));  
        for j = 1:length(un_oris)
            f_ori = find(ori3 == [un_oris(j)]);
            LGN_matallz(i,f_ori,k) = zscore(LGN_mat3(i,f_ori,k));  
            LGN_matshuff(i,f_ori,k) = LGN_matallz(i,f_ori(randperm(length(f_ori))),k);
        end
    end
end
for i = 1:size(LM_mat3,1)
    for k = 1:size(LM_mat3,3)
        LM_matonez(i,:,k) = zscore(LM_mat3(i,:,k));  
        for j = 1:length(un_oris)
            f_ori = find(ori3 == [un_oris(j)]);
            LM_matallz(i,f_ori,k) = zscore(LM_mat3(i,f_ori,k));  
            LM_matshuff(i,f_ori,k) = LM_matallz(i,f_ori(randperm(length(f_ori))),k);
        end
    end
end
%%
load('p9_LGN_spikes_filters.mat')

switch n
    case 1
        LGN_mat = LGN_mat(:,locind,:);
        LGN_mat = LGN_mat(:,locind_match,:);
    case 0
        LGN_mat = LGN_mat(:,awakeind,:);
        LGN_mat = LGN_mat(:,awakeind_match,:);
end

LGN_mat2 = [];
j = 1;
for i = 1:winstep:(size(LGN_mat,3)-winwidth)
    LGN_mat2(:,:,j) = sum(LGN_mat(:,:,i:i+(winwidth-1)),3);
    j = j+1;
end

f_ori = [];
for j = 1:length(un_oris)
    f_ori = [f_ori find(ori2 == [un_oris(j)])];
end
f_ori = sort(f_ori);
ori3 = ori2(f_ori);
LGN_mat3 = LGN_mat2(:,f_ori,:);

LGN_matshuff = zeros(size(LGN_mat3,1),size(LGN_mat3,2),size(LGN_mat3,3));
LGN_matonez = zeros(size(LGN_mat3,1),size(LGN_mat3,2),size(LGN_mat3,3));
LGN_matallz = zeros(size(LGN_mat3,1),size(LGN_mat3,2),size(LGN_mat3,3));
for i = 1:size(LGN_mat3,1)
    for k = 1:size(LGN_mat3,3)
        LGN_matonez(i,:,k) = zscore(LGN_mat3(i,:,k));  
        for j = 1:length(un_oris)
            f_ori = find(ori3 == [un_oris(j)]);
            LGN_matallz(i,f_ori,k) = zscore(LGN_mat3(i,f_ori,k));  
            LGN_matshuff(i,f_ori,k) = LGN_matallz(i,f_ori(randperm(length(f_ori))),k);
        end
    end
end

%%
load('p9_LM_spikes_filters.mat')

switch n
    case 1
        LM_mat = LM_mat(:,locind,:);
        LM_mat = LM_mat(:,locind_match,:);
    case 0
        LM_mat = LM_mat(:,awakeind,:);
        LM_mat = LM_mat(:,awakeind_match,:);
end

LM_mat2 = [];
j = 1;
for i = 1:winstep:(size(LM_mat,3)-winwidth)
    LM_mat2(:,:,j) = sum(LM_mat(:,:,i:i+(winwidth-1)),3);
    j = j+1;
end

f_ori = [];
for j = 1:length(un_oris)
    f_ori = [f_ori find(ori2 == [un_oris(j)])];
end
f_ori = sort(f_ori);
ori3 = ori2(f_ori);
LM_mat3 = LM_mat2(:,f_ori,:);

LM_matshuff = zeros(size(LM_mat3,1),size(LM_mat3,2),size(LM_mat3,3));
LM_matonez = zeros(size(LM_mat3,1),size(LM_mat3,2),size(LM_mat3,3));
LM_matallz = zeros(size(LM_mat3,1),size(LM_mat3,2),size(LM_mat3,3));
for i = 1:size(LM_mat3,1)
    for k = 1:size(LM_mat3,3)
        LM_matonez(i,:,k) = zscore(LM_mat3(i,:,k));  
        for j = 1:length(un_oris)
            f_ori = find(ori3 == [un_oris(j)]);
            LM_matallz(i,f_ori,k) = zscore(LM_mat3(i,f_ori,k));  
            LM_matshuff(i,f_ori,k) = LM_matallz(i,f_ori(randperm(length(f_ori))),k);
%         LM_mat3(i,:,k) = zscore(LM_mat3(i,:,k));  
%         LM_matshuff(i,:,k) = LM_mat3(i,randperm(size(LM_mat3,2)),k);
        end
    end
end


%%
if exist('c','var') == 0
    c = cvpartition(length(ori3),'KFold',10);
end
A = zeros(10,size(LGN_mat3,1),nsub,size(LGN_mat3,3),(200/winstep)+1);
B = zeros(10,size(V_mat3,1),nsub,size(V_mat3,3),(200/winstep)+1);
% UV = zeros(10,length(ori3),10,size(V_mat3,3),201);
r = [];
for h = 1:10
    for i = 1:size(LGN_mat3,3)
        for j = (-100/winstep):(100/winstep)
            if i+j > 0 && i+j < size(V_mat3,3)
                [Atemp,Btemp,rtemp,Utemp,Vtemp] = canoncorr(transpose(squeeze(LGN_matshuff(:,training(c,h),i))),transpose(squeeze(V_matshuff(:,training(c,h),i+j))));
                A(h,:,:,i,j+21) = Atemp(:,1:nsub);
                B(h,:,:,i,j+21) = Btemp(:,1:nsub);
                testsize = sum(test(c,h));
                trainsize = sum(training(c,h));

                Utrain = transpose(squeeze(A(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(LGN_matonez(:,training(c,h),i));
                Vtrain = transpose(squeeze(B(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,training(c,h),i+j));
                UVtrain = transpose(vertcat(Utrain,Vtrain));
                Utest = transpose(squeeze(A(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(LGN_matonez(:,test(c,h),i));
                Vtest = transpose(squeeze(B(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,test(c,h),i+j));
                UVtest = transpose(vertcat(Utest,Vtest));
                cc = corrcoef(Utest(1,:),Vtest(1,:));
                cc(1,2);
                r(h,i,j+(100/winstep)+1) = cc(1,2);
                Umdllda = fitcdiscr(transpose(Utrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Umdllda, transpose(Utest));
                pu(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                Vmdllda = fitcdiscr(transpose(Vtrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Vmdllda, transpose(Vtest));
                pv(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;   

                UVmdllda = fitcdiscr(UVtrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(UVmdllda, UVtest);
                p(h,i,j+21) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;  
                p(h,i,j+21)

            else 
                r(h,i,j+(100/winstep)+1) = NaN;
                p(h,i,j+(100/winstep)+1) = NaN;
                pu(h,i,j+(100/winstep)+1) = NaN;
                pv(h,i,j+(100/winstep)+1) = NaN;

            end
        end
    end
end
ravg = squeeze(mean(r,1,'omitnan'));
ravg(ravg < 0) = NaN;
pavg = squeeze(mean(p,1));
% pavg2 = squeeze(sum(p2,1)/length(ori3));
%%
figure;
h = heatmap(squeeze(ravg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((size(V_mat,3)-winwidth-50):-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
h.Title = 'LGN-V1 CCA Correlation Coefficiencts (Stationary, orientations calculated together';
switch n
    case 0
        h.Title = 'LGN-V1 CCA Correlation Coefficients (Stationary, orientations calculated together';
        savefig(fullfile(pwd,'figs',strcat('LGNV1_CCAcorrelations_stat_together_',ext)));
    case 1
        h.Title = 'LGN-V1 CCA Correlation Coefficients (Locomoting, orientations calculated together';
        savefig(fullfile(pwd,'figs',strcat('LGNV1_CCAcorrelations_loc_together_',ext)));
end
figure;
h = heatmap(squeeze(pavg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((size(V_mat,3)-winwidth-50):-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
switch n
    case 0
        h.Title = 'LGN-V1 CCA Decoding (Stationary, orientations calculated together';
        savefig(fullfile(pwd,'figs',strcat('LGNV1_CCAdecoding_stat_together_',ext)));
    case 1
        h.Title = 'LGN-V1 CCA Decoding (Locomoting, orientations calculated together';
        savefig(fullfile(pwd,'figs',strcat('LGNV1_CCAdecoding_loc_together_',ext)));
end
switch n
    case 0
        save(strcat(animal,'_LGNV1_CCA_stat_together_',ext,'.mat'),'c','LGN_matonez','V_matonez','LGN_matallz','V_matallz','ori3','r','p','pu','pv','A','B','awakeind','awakeind_match','-v7.3');
    case 1
        save(strcat(animal,'_LGNV1_CCA_loc_together_',ext,'.mat'),'c','LGN_matonez','V_matonez','LGN_matallz','V_matallz','ori3','r','p','pu','pv','A','B','locind','locind_match','-v7.3');
end

% figure;
% h = heatmap(squeeze(pavg2(:,:)),'Colormap',parula,'GridVisible','off');
% h.YDisplayData = flipud(h.YDisplayData);
% h.XDisplayLabels = nan(size(h.XDisplayData));
% h.YDisplayLabels = nan(size(h.YDisplayData));
% h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
% h.YDisplayLabels(1:(50/winstep):size(h.YDisplayLabels)) = num2cell(800:-(50):0);
% h.XLabel = 'Delay';
% h.YLabel = 'Time from Stimulus Onset';
% switch n
%     case 0
%         h.Title = 'LGN-V1 Full Neural Activity Decoding Heatmap (Stationary)';
%         savefig(fullfile(pwd,'figs/LGNV1_stat_fullpop'));
%     case 1
%         h.Title = 'LGN-V1 Full Neural Activity Decoding Heatmap (Locomoting)';
%         savefig(fullfile(pwd,'figs/LGNV1_loc_fullpop'));
% end
% switch n
%     case 0
%         save(strcat(animal,'_LGNV1_stat_fullpop.mat'),'c','p2','V_matonez','LGN_matonez','awakeind','awakeind_match','-v7.3');
%     case 1
%         save(strcat(animal,'_LGNV1_loc_fullpop.mat'),'c','p2','V_matonez','LGN_matonez','locind','locind_match','-v7.3');
% end
%%
A = zeros(10,size(V_mat3,1),nsub,size(V_mat3,3),(200/winstep)+1);
B = zeros(10,size(LM_mat3,1),nsub,size(LM_mat3,3),(200/winstep)+1);
r = [];
for h = 1:10
    for i = 1:size(V_mat3,3)
        for j = (-100/winstep):(100/winstep)
            if i+j > 0 && i+j < size(LM_mat3,3)
                testsize = sum(test(c,h));
                [Atemp,Btemp,rtemp,Utemp,Vtemp] = canoncorr(transpose(squeeze(V_matshuff(:,training(c,h),i))),transpose(squeeze(LM_matshuff(:,training(c,h),i+j))));
                A(h,:,:,i,j+21) = Atemp(:,1:nsub);
                B(h,:,:,i,j+21) = Btemp(:,1:nsub);
                Utrain = transpose(squeeze(A(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,training(c,h),i));
                Vtrain = transpose(squeeze(B(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(LM_matonez(:,training(c,h),i+j));
                UVtrain = transpose(vertcat(Utrain,Vtrain));
                Utest = transpose(squeeze(A(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(V_matonez(:,test(c,h),i));
                Vtest = transpose(squeeze(B(h,:,1:nsub,i,j+(100/winstep)+1))) * squeeze(LM_matonez(:,test(c,h),i+j));
                UVtest = transpose(vertcat(Utest,Vtest));
                cc = corrcoef(Utest(1,:),Vtest(1,:));
                cc(1,2);
                r(h,i,j+(100/winstep)+1) = cc(1,2);

                Umdllda = fitcdiscr(transpose(Utrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Umdllda, transpose(Utest));
                pu(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;

                Vmdllda = fitcdiscr(transpose(Vtrain),ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(Vmdllda, transpose(Vtest));
                pv(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;   

                UVmdllda = fitcdiscr(UVtrain,ori3(training(c,h)),'DiscrimType','pseudolinear');
                ptemp = predict(UVmdllda, UVtest);
                p(h,i,j+(100/winstep)+1) = sum(transpose(ptemp) == ori3(test(c,h)))/testsize;  
                p(h,i,j+(100/winstep)+1)
            else 
                r(h,i,j+(100/winstep)+1) = NaN;
                p(h,i,j+(100/winstep)+1) = NaN;
                pu(h,i,j+(100/winstep)+1) = NaN;
                pv(h,i,j+(100/winstep)+1) = NaN;

            end
        end
    end
end
ravg = squeeze(mean(r,1,'omitnan'));
ravg(ravg < 0) = NaN;
pavg = squeeze(mean(p,1));
%%
figure;
h = heatmap(squeeze(ravg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((1000-winwidth-50):-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
switch n
    case 0
        h.Title = 'V1-LM CCA Correlation Coefficiencts (Stationary, orientations calculated together';
        savefig(fullfile(pwd,'figs',strcat('V1LM_CCAcorrelations_stat_together_',ext)));
    case 1
        h.Title = 'V1-LM CCA Correlation Coefficiencts (Locomoting, orientations calculated together';
        savefig(fullfile(pwd,'figs',strcat('V1LM_CCAcorrelations_loc_together_',ext)));
end

figure;
h = heatmap(squeeze(pavg(:,:)),'Colormap',parula,'GridVisible','off');
h.YDisplayData = flipud(h.YDisplayData);
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
h.YDisplayLabels(5:(50/winstep):size(h.YDisplayLabels)) = num2cell((size(V_mat,3)-winwidth-50):-(50):0);
h.XLabel = 'Delay';
h.YLabel = 'Time from Stimulus Onset';
switch n
    case 0
        h.Title = 'V1-LM CCA Decoding (Stationary, orientations calculated together';
        savefig(fullfile(pwd,'figs',strcat('LGNV1_CCAdecoding_stat_together_',ext)));
    case 1
        h.Title = 'V1-LM CCA Decoding (Locomoting, orientations calculated together';
        savefig(fullfile(pwd,'figs',strcat('V1LM_CCAdecoding_loc_together_',ext)));
end
switch n
    case 0
        save(strcat(animal,'_V1LM_CCA_stat_together_',ext,'.mat'),'c','V_matonez','LM_matonez','V_matallz','LM_matallz','ori3','r','p','pu','pv','A','B','awakeind','awakeind_match','-v7.3');
    case 1
        save(strcat(animal,'_V1LM_CCA_loc_together_',ext,'.mat'),'c','V_matonez','LM_matonez','V_matallz','LM_matallz','ori3','r','p','pu','pv','A','B','locind','locind_match','-v7.3');
end

% figure;
% h = heatmap(squeeze(pavg2(:,:)),'Colormap',parula,'GridVisible','off');
% h.YDisplayData = flipud(h.YDisplayData);
% h.XDisplayLabels = nan(size(h.XDisplayData));
% h.YDisplayLabels = nan(size(h.YDisplayData));
% h.XDisplayLabels(1:(size(h.XDisplayLabels)-1)/4:size(h.XDisplayLabels)) = num2cell(-100:50:100);
% h.YDisplayLabels(1:(50/winstep):size(h.YDisplayLabels)) = num2cell(800:-(50):0);
% h.XLabel = 'Delay';
% h.YLabel = 'Time from Stimulus Onset';
% h.Title = 'V1-LM Full Neural Activity Decoding Heatmap (Locomoting)';
% switch n
%     case 0
%         h.Title = 'V1-LM Full Neural Activity Decoding Heatmap (Stationary)';
%         savefig(fullfile(pwd,'figs/V1LM_stat_fullpop'));
%     case 1
%         h.Title = 'V1-LM Full Neural Activity Decoding Heatmap (Locomoting)';
%         savefig(fullfile(pwd,'figs/V1LM_loc_fullpop'));
% end
% switch n
%     case 0
%         save(strcat(animal,'_V1LM_stat_fullpop.mat'),'c','p2','V_matonez','LM_matonez','awakeind','awakeind_match','-v7.3');
%     case 1
%         save(strcat(animal,'_V1LM_loc_fullpop.mat'),'c','p2','V_matonez','LM_matonez','locind','locind_match','-v7.3');
% end
% scatter(U(:,1),V(:,1))



%% A(,:1) * test = test data mapped onto first canonical variable. Repeat for all test trials. Correlate test trials between areas. 