agl_pos = -3;
step = 3;
fit_record = NaN*ones(3,8);
%1st col: angular position
%2nd col: negative angular velocity (for 2 Gaussian fit)
%3rd col: positive angular velocity (for 2 Gaussian fit)
%otherwise 2nd col = 3rd col 
%4th col: standard deviation of 2nd col
%5th col: standard deviation of 3rd col
%6th col: proportion of Gaussian of negative angular velocity
%7th col: proportion of Gaussian of positive angular velocity 
%8th col: effective angular velocity 

torque_data = cell(3,2);
%for torque analysis section 
%1st col: angular position
%2nd col: statistical data based on histogram, calculated based on
%fluctuation theorem 
%inside 2nd col: 1st col: delta 2theta/kT

hist_record = cell(3,3);
%to store histogram data of angular velocity
%1st col: angular position
%2nd col: angular velocity bin values
%3rd col: counts of each angular velocity bin

numGauss = cell(3,3);
%store results of fitting number of Gaussian components based on AIC and BIC 
%1st col: angular position
%2nd col: 1 Gaussian
%3rd col: 2 Gaussian


i = 1;
while agl_pos <= 130
    
    [numComponents, MU, COV, PPp, histdata] = GMM(d_ang_raw, agl_pos);  
    %[numComponents, MU, COV, PPp, histdata] = GMM(d_ang_raw_dwell_shift, agl_pos);  
    %for 42 and 45, need to remove rng('default') and rng(1) (for d_ang_raw_dwell_shift)
    
    %for torque, use this 
    %[numComponents, MU, COV, PPp, torque_FT] = GMM2(d_ang_raw, agl_pos);      
    %[numComponents, MU, COV, PPp, torque_FT] = GMM2(d_ang_raw_dwell_shift, agl_pos);  
    
    fit_record(i,1) = agl_pos;    
    fit_record(i,2:3) = MU;
    fit_record(i,4:5) = sqrt(COV); 
    fit_record(i,6:7) = PPp; 
    %torque_data{i,1} = agl_pos;
    %torque_data{i,2} = torque_FT; 
    hist_record{i,1} = agl_pos;
    hist_record(i,2:3) = histdata(1:2,1);
    %hist_record{i,2} = histdata{1,1};
    %hist_record{i,3} = histdata{2,1};
    
    numGauss{i,1} = agl_pos;
    numGauss(i,2:3)= histdata(3:4,1);
    
    if numComponents == 2
        fit_record(i,8) = sum(fit_record(i,2:3).*fit_record(i,6:7));
    else 
        fit_record(i,8) = fit_record(i,2);
    end
    %check to make sure 2nd col stores only negative angular velocity
    %and other col follows the switch 
    if fit_record(i,3) < fit_record(i,2)
        temp = fit_record(i,2) ;
        fit_record(i,2) = fit_record(i,3);
        fit_record(i,3) = temp;
        
        temp2 = fit_record(i,4);
        fit_record(i,4) = fit_record(i,5);
        fit_record(i,5) = temp2;
        
        temp3 = fit_record(i,6);
        fit_record(i,6) = fit_record(i,7);
        fit_record(i,7) = temp3;
    end
    
    i = i+1;
    agl_pos = agl_pos + step;
end
close all

%plot angular velocity
figure('pos',[1 1 948 605])
plot(fit_record(:,1),fit_record(:,2),'b.-',fit_record(:,1),fit_record(:,3),'r.-')
hold on
plot(fit_record(:,1),fit_record(:,8),'k.-')
legend('Negative angular velocity','Positive angular velocity','Average angular velocity')
xlabel('Angular position (deg)')
ylabel('Angular velocity (deg/ms)');
title('Average angular velocity - angular position')
xlim([-4 131])

%Proportion of each components of double Gaussian
temp = find(fit_record(:,6) ~= 1);

figure
subplot(2,1,1)
bar(fit_record(temp,1),fit_record(temp,6:7),'stacked');
legend('Negative','Positive');
xlabel('Angular position (deg)')
ylabel('Proportion');
title('Proportion distribution')
xlim([35 70])

subplot(2,1,2)
plot(fit_record(temp,1),fit_record(temp,7) - fit_record(temp,6),'.')
xlabel('Angular position (deg)')
ylabel('Proportion difference');
title('Difference b/w positive and negative proportion')

%Standard deviation of 2 Gaussian 
figure
h_bw = scatter(fit_record(temp,1),fit_record(temp,4),'b')
h_1 = h_bw;
hold on
h_fw = scatter(fit_record(temp,1),fit_record(temp,5),'r');

%standard deviation of fitting at remaining angle of 1 Gaussian
h_G = scatter(fit_record(fit_record(:,6) == 1,1), fit_record(fit_record(:,6) == 1,4),'k');

[h, ~] = legend([h_1, h_fw, h_G],{'Backward','Forward','1 Gaussian'});

xlabel('Angle (deg)');
ylabel('Standard deviation (deg/ms)')
title('Standard deviation of angular velocity distribution');
xlim([-5 130]);

%% Plot 1 vs 2 Gaussian 

for i = 1 : length(numGauss)
    %for laptop
    %fpath = 'D:\Onedrive-NTU\OneDrive - Nanyang Technological University\PhD\rev analyze\Results\Ang vel';
    %for desktop at school 
    fpath = 'D:\Onedrive-Home\OneDrive - Nanyang Technological University\PhD\rev analyze\Results\Ang vel'; 
    
    avel = hist_record{i,2}; cvel = hist_record{i,3};
    agl_pos = numGauss{i,1};
    g = figure
    for numComponents = 1 : 2              
        
        switch numComponents
            case 1
                 paramEsts = numGauss{i,2};
                 MU = paramEsts.mu;
                 COV = paramEsts.Sigma;   %sigma here is the covariance 
                 PPp = paramEsts.PComponents;
            case 2
                 paramEsts = numGauss{i,3};
                MU=[paramEsts.mu(1);paramEsts.mu(2);];
                COV = cat(3,[paramEsts.Sigma(1)],[paramEsts.Sigma(2)]);  %store covariance 
                PPp = [paramEsts.PComponents(1),paramEsts.PComponents(2)];  %component proportion
        end
        
        objA = gmdistribution(MU,COV,PPp);
        xgridss=transpose(linspace(min(avel),max(avel),100)); 
        
        subplot(2,1,numComponents);
        
        plot(avel,cvel,'ko'); hold on 
        plot(xgridss,pdf(objA,xgridss)*sum(cvel*binwidth),'r-','linewidth',2)  %convert pdf to real statistic    
        if numComponents == 2  
            obj_1 = gmdistribution(MU(1), COV(1), PPp(1));
            obj_2 = gmdistribution(MU(2), COV(2), PPp(2));
            plot(xgridss,pdf(obj_1,xgridss)*sum(cvel*binwidth)*PPp(1))
            plot(xgridss,pdf(obj_2,xgridss)*sum(cvel*binwidth)*PPp(2))
            %plot(xgridss,pdf(obj_1,xgridss)*sum(cvel*binwidth)*PPp(1) + pdf(obj_2,xgridss)*sum(cvel*binwidth)*PPp(2))
        end
        title(sprintf('Angular velocity distribution at %d deg - %d Gaussian',agl_pos, numComponents))
        xlabel('Angular velocity (deg/ms)')
        ylabel('Counts')      
    end
    saveas(gca, fullfile(fpath, sprintf('angle = %d',agl_pos)), 'png');   
    saveas(gca, fullfile(fpath, sprintf('angle = %d',agl_pos)), 'fig'); 
    
    
    
    
    
end
close all

%% Smooth angular velocity - angular position curve
data_need = {'fit_record','fit_record_lim','fit_record_rawdwshift_lim','fit_record_rawdwshift'};

for i = 1 : length(data_need)
    load(data_need{i});
end

%plot all 4 together 
figure
plot(fit_record(:,1),fit_record(:,2),'b.-',fit_record(:,1),fit_record(:,3),'r.-')
hold on
plot(fit_record(:,1),fit_record(:,8),'k.-')

plot(fit_record_lim(:,1),fit_record_lim(:,2),'b.-',fit_record_lim(:,1),fit_record_lim(:,3),'r.-')
hold on
plot(fit_record_lim(:,1),fit_record_lim(:,8),'k.-')

plot(fit_record_rawdwshift(:,1),fit_record_rawdwshift(:,2),'b.-',fit_record_rawdwshift(:,1),fit_record_rawdwshift(:,3),'r.-')
hold on
plot(fit_record_rawdwshift(:,1),fit_record_rawdwshift(:,8),'k.-')

plot(fit_record_rawdwshift_lim(:,1),fit_record_rawdwshift_lim(:,2),'b.-',fit_record_rawdwshift_lim(:,1),fit_record_rawdwshift_lim(:,3),'r.-')
hold on
plot(fit_record_rawdwshift_lim(:,1),fit_record_rawdwshift_lim(:,8),'k.-')

xlabel('Angular position (deg)')
ylabel('Angular velocity (deg/ms)');
title('Average angular velocity - angular position')
xlim([-4 121])
%legend('No lim','Lim','Shift-no Lim','Shift-Lim')
%obtain average and plot
angvel_ave = NaN*ones(length(fit_record),4);
angvel_ave(:,1) = fit_record(:,1);

bin = 1.645^2;
for i = 1:length(fit_record)
   %mean
   angvel_ave(i,2) = (fit_record(i,2) + fit_record_lim(i,2) + fit_record_rawdwshift(i,2) + fit_record_rawdwshift_lim(i,2))/4;
   angvel_ave(i,3) = (fit_record(i,3) + fit_record_lim(i,3) + fit_record_rawdwshift(i,3) + fit_record_rawdwshift_lim(i,3))/4;
   angvel_ave(i,4) = (fit_record(i,8) + fit_record_lim(i,8) + fit_record_rawdwshift(i,8) + fit_record_rawdwshift_lim(i,8))/4;
   %standard deviation of 1 Gaussian 
   %based on https://stats.stackexchange.com/questions/241678/confidence-interval-on-sum-of-estimates-vs-estimate-of-whole?noredirect=1&lq=1
   % (Ex 4) https://www.statlect.com/probability-distributions/normal-distribution-linear-combinations
   
   angvel_ave(i,5) = sqrt(fit_record(i,4)^2 + fit_record_lim(i,4)^2 + fit_record_rawdwshift(i,4)^2 + fit_record_rawdwshift_lim(i,4)^2)/4;
   angvel_ave(i,6) = sqrt(fit_record(i,5)^2 + fit_record_lim(i,5)^2 + fit_record_rawdwshift(i,5)^2 + fit_record_rawdwshift_lim(i,5)^2)/4;
   %90% confidence interval of 1 Gaussian (z = 1.645)
   angvel_ave(i,5) = 1.645*angvel_ave(i,5)/sqrt(bin);
   angvel_ave(i,6) = 1.645*angvel_ave(i,6)/sqrt(bin);
   %calculate mean standard deviation of 2 Gaussian for each data set 
   fit_record(i,8) = sqrt(fit_record(i,6)^2*fit_record(i,4)^2 + (1-fit_record(i,6))^2*fit_record(i,5)^2);
   fit_record_lim(i,8) = sqrt(fit_record_lim(i,6)^2*fit_record_lim(i,4)^2 + (1-fit_record_lim(i,6))^2*fit_record_lim(i,5)^2);
   fit_record_rawdwshift(i,8) = sqrt(fit_record_rawdwshift(i,6)^2*fit_record_rawdwshift(i,4)^2 + (1-fit_record_rawdwshift(i,6))^2*fit_record_rawdwshift(i,5)^2);
   fit_record_rawdwshift_lim(i,8) = sqrt(fit_record_rawdwshift_lim(i,6)^2*fit_record_rawdwshift_lim(i,4)^2 + (1-fit_record_rawdwshift_lim(i,6))^2*fit_record_rawdwshift_lim(i,5)^2);
   %calcuate mean standard deviation
   angvel_ave(i,7) = sqrt(fit_record(i,8)^2 + fit_record_lim(i,8)^2 + fit_record_rawdwshift(i,8)^2 + fit_record_rawdwshift_lim(i,8)^2)/4;
   
   %confidence interval of 1 Gaussian 
   angvel_ave(i,7) = 1.645*angvel_ave(i,7)/sqrt(bin);
       
end

%plot all 4 together 
temp = find(fit_record(:,6) ~= 1);

figure
errorbar(angvel_ave(:,1),angvel_ave(:,4),angvel_ave(:,7),'k')
hold on
errorbar(angvel_ave(temp,1),angvel_ave(temp,3),angvel_ave(temp,6),'r')
errorbar(angvel_ave(temp,1),angvel_ave(temp,2),angvel_ave(temp,5),'b')
xlim([-4 121])
legend('Average','Positive','Negative')
xlabel('Angular position (deg)')
ylabel('Angular velocity (deg/ms)');
title('Average angular velocity - angular position')

figure
plot(angvel_ave(temp,1),angvel_ave(temp,2),'b.-',angvel_ave(temp,1),angvel_ave(temp,3),'r.-')
hold on
plot(angvel_ave(:,1),angvel_ave(:,4),'k.-')

%smooth data 
%sm_rloess = smooth(angvel_ave(:,1),angvel_ave(:,4),0.25,'rloess');
%sm_rlowess = smooth(angvel_ave(:,1),angvel_ave(:,4),0.1,'rlowess');
sm_sgolay= smooth(angvel_ave(:,1),angvel_ave(:,4),0.1,'sgolay');
%plot(angvel_ave(:,1),angvel_ave(:,4),'k.-')
%hold on
%plot(angvel_ave(:,1),sm_rloess,'c.-')
%plot(angvel_ave(:,1),sm_rlowess,'b.-')
%plot(angvel_ave(:,1),sm_sgolay,'r.-')

sm2_sgolay= smooth(angvel_ave(:,1),angvel_ave(:,2),0.1,'sgolay');
sm3_sgolay= smooth(angvel_ave(:,1),angvel_ave(:,3),0.1,'sgolay');

plot(angvel_ave(:,1),angvel_ave(:,4),'k.-')
hold on
plot(angvel_ave(temp,1),sm2_sgolay(temp),'.-','Color',[0.2 .4 0.5])
plot(angvel_ave(temp,1),sm3_sgolay(temp),'.-','Color', [0.8 0.7 0.3])
plot(angvel_ave(:,1),sm_sgolay,'r.-')
xlabel('Angular position (deg)')
ylabel('Angular velocity (deg/ms)');
title('Average angular velocity - angular position')
xlim([-4 121])
legend('Average from 4 cases','Smooth negative part','Smooth positive part','Smooth overall')

figure
errorbar(angvel_ave(:,1),sm_sgolay,angvel_ave(:,7),'k.-')
hold on
errorbar(angvel_ave(temp,1),sm3_sgolay(temp),angvel_ave(temp,6),'r.-')
errorbar(angvel_ave(temp,1),sm2_sgolay(temp),angvel_ave(temp,5),'b.-')
xlim([-1 121])
legend('Average','Positive','Negative')
xlabel('Angular position (deg)')
ylabel('Angular velocity (deg/ms)');
title('Average angular velocity - angular position')

%follow by torque 
kappa = 50; %pN.nm/rad^2
temp = find(fit_record(:,1) >= 78 & fit_record(:,1) <= 125);
angle = fit_record(temp,1);
aglvel = sm_sgolay(temp);
fit_result = fit(angle,aglvel,'poly1');
plot(fit_result,angle, aglvel);
coeff = coeffvalues(fit_result);
friction = -kappa/coeff(1); %pN.nm.ms;

torque_smooth = friction*sm_sgolay*pi/180;
torque_pos = friction*sm3_sgolay*pi/180;
torque_neg = friction*sm2_sgolay*pi/180;
temp = find(fit_record(:,6) ~= 1);

plot(fit_record(:,1),torque_smooth ,'.-');
hold on
plot(fit_record(temp,1),torque_pos(temp) ,'.-');
plot(fit_record(temp,1),torque_neg(temp) ,'.-');

figure
errorbar(angvel_ave(:,1),torque_smooth,friction*angvel_ave(:,7)*pi/180,'k')
hold on
errorbar(angvel_ave(temp,1),torque_pos(temp),friction*angvel_ave(temp,6)*pi/180,'r')
errorbar(angvel_ave(temp,1),torque_neg(temp),friction*angvel_ave(temp,5)*pi/180,'b')
xlim([-1 121])
legend('Average','Positive','Negative')
xlabel('Angular position (deg)')
ylabel('Torque (pN.nm)');
title('Average torque - angular position')

%calculate the energy 
work = 0;
for i = 1 : length(torque_smooth)-1
    work = work + sum(torque_smooth(i:i+1))*3*pi/(2*180);
end

%% Smooth torque - position profile
%torque from Fluctuation theorem 
data_need = {'fit_record','fit_record_lim','fit_record_rawdwshift_lim','fit_record_rawdwshift'};

for i = 1 : length(data_need)
    load(data_need{i});
end
R = 20; %nm
rc = 19; %nm
eta = 1e-9; %pN.s/nm^2
friction = 8*pi*eta*R^3 + 6*pi*eta*R*rc^2;

torque_profile = NaN*ones(length(fit_record),7);
torque_profile(:,1) = fit_record(:,1);
temp = find(fit_record(:,6) ~= 1);
%1st col: angle
%2nd col: average torque calculated from fit_record
%3rd col: average torque calculated from fit_record_lim
%4th col: average torque calculated from fit_record_rawdwshift
%5th col: average torque calculated from fit_record_rawdwshift_lim
%6th col: average of 2nd - 5th col
%7th col: std of average torque of 6th col
torque_pos = NaN*ones(length(fit_record),7);
%positive torque
%convention like torque_profile, but for positive torque
torque_pos(:,1) = fit_record(:,1);
torque_neg = NaN*ones(length(fit_record),7);
%negative torque
%convention like torque_profile, but for negative torque
torque_neg(:,1) = fit_record(:,1);

for i = 1 : length(fit_record)
    %Convert angular velocity and standard deviation from deg/ms to rad 
    torque_temp1 = NaN*ones(1,5);
    torque_temp2 = NaN*ones(1,5);
    torque_temp3 = NaN*ones(1,5);
    torque_temp4 = NaN*ones(1,5);
    %1st col: neg torque
    %2nd col: pos torque
    %3rd col: std of neg torque
    %4th col: std of pos torque
    %5th col: std of average torque *1000 to convert from rad/ms to rad/s 
    %Convert angular velocity and standard deviation from deg/ms to rad 

    torque_temp1(1,1) = 2*4.14*(fit_record(i,2)*pi/(180*100))/(fit_record(i,4)*pi/(180*100))^2;   %1kT = 4.14 pN.nm
    torque_temp1(1,2) = 2*4.14*(fit_record(i,3)*pi/(180*100))/(fit_record(i,5)*pi/(180*100))^2;    %1kT = 4.14 pN.nm
    torque_temp1(1,3) = friction*fit_record(i,4)*pi*1000/180; %error of neg torque 
    torque_temp1(1,4) = friction*fit_record(i,5)*pi*1000/180; %error of pos torque 

    
    torque_temp2(1,1) = 2*4.14*(fit_record_lim(i,2)*pi/(180*100))/(fit_record_lim(i,4)*pi/(180*100))^2;   %1kT = 4.14 pN.nm
    torque_temp2(1,2) = 2*4.14*(fit_record_lim(i,3)*pi/(180*100))/(fit_record_lim(i,5)*pi/(180*100))^2;    %1kT = 4.14 pN.nm
    torque_temp2(1,3) = friction*fit_record_lim(i,4)*pi*1000/180;
    torque_temp2(1,4) = friction*fit_record_lim(i,5)*pi*1000/180;

    torque_temp3(1,1) = 2*4.14*(fit_record_rawdwshift(i,2)*pi/(180*100))/(fit_record_rawdwshift(i,4)*pi/(180*100))^2;   %1kT = 4.14 pN.nm
    torque_temp3(1,2) = 2*4.14*(fit_record_rawdwshift(i,3)*pi/(180*100))/(fit_record_rawdwshift(i,5)*pi/(180*100))^2;    %1kT = 4.14 pN.nm
    torque_temp3(1,3) = friction*fit_record_rawdwshift(i,4)*pi*1000/180;
    torque_temp3(1,4) = friction*fit_record_rawdwshift(i,5)*pi*1000/180;

    
    torque_temp4(1,1) = 2*4.14*(fit_record_rawdwshift_lim(i,2)*pi/(180*100))/(fit_record_rawdwshift_lim(i,4)*pi/(180*100))^2;   %1kT = 4.14 pN.nm
    torque_temp4(1,2) = 2*4.14*(fit_record_rawdwshift_lim(i,3)*pi/(180*100))/(fit_record_rawdwshift_lim(i,5)*pi/(180*100))^2;    %1kT = 4.14 pN.nm
    torque_temp4(1,3) = friction*fit_record_rawdwshift_lim(i,4)*pi*1000/180;
    torque_temp4(1,4) = friction*fit_record_rawdwshift_lim(i,5)*pi*1000/180;

    torque_pos(i,2) = torque_temp1(1,2); %+ torque(i,3)*fit_record(i,7) 
    torque_pos(i,3) = torque_temp2(1,2);
    torque_pos(i,4) = torque_temp3(1,2);
    torque_pos(i,5) = torque_temp4(1,2);
    torque_pos(i,6) = mean(torque_pos(i,2:5));
    %torque_pos(i,7) = std(torque_pos(i,2:5)); 
    torque_pos(i,7) = sqrt(torque_temp1(1,3)^2 + torque_temp2(1,3)^2 + torque_temp3(1,3)^2 + torque_temp4(1,3)^2)/4;

    torque_neg(i,2) = torque_temp1(1,1); %+ torque(i,3)*fit_record(i,7) 
    torque_neg(i,3) = torque_temp2(1,1);
    torque_neg(i,4) = torque_temp3(1,1);
    torque_neg(i,5) = torque_temp4(1,1);
    torque_neg(i,6) = mean(torque_neg(i,2:5));    
    %torque_neg(i,7) = std(torque_neg(i,2:5));
    torque_neg(i,7) = sqrt(torque_temp1(1,4)^2 + torque_temp2(1,4)^2 + torque_temp3(1,4)^2 + torque_temp4(1,4)^2)/4;

    
    if fit_record(i,6) == fit_record(i,7)
        torque_profile(i,2) = torque_temp1(1,1); %+ torque(i,3)*fit_record(i,7) 
        torque_profile(i,3) = torque_temp2(1,1);
        torque_profile(i,4) = torque_temp3(1,1);
        torque_profile(i,5) = torque_temp4(1,1);
        
        torque_temp1(1,5) = sqrt((fit_record(i,6)*torque_temp1(1,3))^2); %error of average torque 
        torque_temp2(1,5) = sqrt((fit_record_lim(i,6)*torque_temp1(1,3))^2); %error of average torque 
        torque_temp3(1,5) = sqrt((fit_record_rawdwshift(i,6)*torque_temp1(1,3))^2); %error of average torque 
        torque_temp4(1,5) = sqrt((fit_record_rawdwshift_lim(i,6)*torque_temp1(1,3))^2); %error of average torque 
        
        torque_profile(i,7) = sqrt(torque_temp1(1,5)^2 + torque_temp2(1,5)^2 + torque_temp3(1,5)^2 + torque_temp4(1,5)^2)/4;
    else
        torque_profile(i,2) = sum(torque_temp1(1,1:2).*fit_record(i,6:7)); %+ torque(i,3)*fit_record(i,7)    
        torque_profile(i,3) = sum(torque_temp2(1,1:2).*fit_record_lim(i,6:7)); %+ torque(i,3)*fit_record(i,7)    
        torque_profile(i,4) = sum(torque_temp3(1,1:2).*fit_record_rawdwshift(i,6:7)); %+ torque(i,3)*fit_record(i,7)    
        torque_profile(i,5) = sum(torque_temp4(1,1:2).*fit_record_rawdwshift_lim(i,6:7)); %+ torque(i,3)*fit_record(i,7)    
        
        torque_temp1(1,5) = sqrt((fit_record(i,6)*torque_temp1(1,3))^2+(fit_record(i,7)*torque_temp1(1,4))^2); %error of average torque 
        torque_temp2(1,5) = sqrt((fit_record_lim(i,6)*torque_temp2(1,3))^2+(fit_record_lim(i,7)*torque_temp2(1,4))^2); %error of average torque 
        torque_temp3(1,5) = sqrt((fit_record_rawdwshift(i,6)*torque_temp3(1,3))^2+(fit_record_rawdwshift(i,7)*torque_temp3(1,4))^2); %error of average torque 
        torque_temp4(1,5) = sqrt((fit_record_rawdwshift_lim(i,6)*torque_temp4(1,3))^2+(fit_record_rawdwshift_lim(i,7)*torque_temp4(1,4))^2); %error of average torque 
        
        torque_profile(i,7) = sqrt(torque_temp1(1,5)^2 + torque_temp2(1,5)^2 + torque_temp3(1,5)^2 + torque_temp4(1,5)^2)/4;
    end
    
    torque_profile(i,6) = mean(torque_profile(i,2:5));
    %torque_profile(i,7) = std(torque_profile(i,2:5));
end




plot(torque_profile(:,1), torque_profile(:,2:5),'k.-')
hold on
plot(torque_pos(:,1), torque_pos(:,2:5),'r.-')
plot(torque_neg(:,1), torque_neg(:,2:5),'b.-')

%plot the average of 4 data sets
figure
plot(torque_profile(:,1), torque_profile(:,6),'k.-')
hold on
plot(torque_pos(:,1), torque_pos(:,6),'r.-')
plot(torque_neg(:,1), torque_neg(:,6),'b.-')

%torque_smooth = smooth(torque_profile(:,1),torque_profile(:,6),0.1,'sgolay');
%torquepos_smooth= smooth(torque_pos(:,1),torque_pos(:,6),0.1,'sgolay');
%torqueneg_smooth= smooth(torque_neg(:,1),torque_neg(:,6),0.1,'sgolay');

torque_smooth = smooth(torque_profile(:,1),torque_profile(:,6),5,'moving');
torquepos_smooth= smooth(torque_pos(:,1),torque_pos(:,6),5,'moving');
torqueneg_smooth= smooth(torque_neg(:,1),torque_neg(:,6),5,'moving');

%torque_smooth = smooth(torque_profile(:,1),torque_profile(:,6),'loess');
%torquepos_smooth= smooth(torque_pos(:,1),torque_pos(:,6),'loess');
%torqueneg_smooth= smooth(torque_neg(:,1),torque_neg(:,6),'loess');

figure
plot(torque_profile(:,1), torque_smooth,'k.-')
hold on
plot(torque_pos(temp,1), torquepos_smooth(temp),'r.-')
plot(torque_neg(temp,1), torqueneg_smooth(temp),'b.-')
xlabel('Angle (deg)')
ylabel('Torque (pN.nm)')
legend('Average','Positive','Negative')


%include error
figure
errorbar(torque_profile(:,1), torque_smooth, torque_profile(:,7))
hold on
errorbar(torque_pos(temp,1), torquepos_smooth(temp), torque_pos(temp,7))
errorbar(torque_neg(temp,1), torqueneg_smooth(temp), torque_neg(temp,7))
xlim([-4 121])
xlabel('Angle (deg)')
ylabel('Torque (pN.nm)')
legend('Average','Positive','Negative')

%calculate the energy 
work = 0;
for i = 1 : length(torque_smooth)-1
    work = work + sum(torque_smooth(i:i+1))*3*pi/(2*180);
end

%figure
%errorbar(torque_profile(:,1), torque_smooth, 1.645*torque_profile(:,7)/4)
%hold on
%errorbar(torque_pos(temp,1), torquepos_smooth(temp), 1.645*torque_pos(temp,7)/4)
%errorbar(torque_neg(temp,1), torqueneg_smooth(temp), 1.645*torque_neg(temp,7)/4)
%xlim([-4 121])
%ylim([-85 220])
%xlabel('Angle (deg)')
%ylabel('Torque (pN.nm)')
%legend('Average','Positive','Negative')

%% Smooth torque-position 2
%by follow definition torque = friction coeff * angular velocity 
data_need = {'fit_record','fit_record_lim','fit_record_rawdwshift_lim','fit_record_rawdwshift'};

for i = 1 : length(data_need)
    load(data_need{i});
end

%value based on
% Nature volume 410, pages 898–904 (19 April 2001)
R = 40; %nm %radius of bead 
rc = 30; %nm %centroid length
eta = 1e-9; %pN.s/nm^2  %water viscosity 
friction = 8*pi*eta*R^3 + 6*pi*eta*R*rc^2;
temp = find(fit_record(:,6) ~= 1);

%calculate average torque with error bar

torque_profile = NaN*ones(length(fit_record),6);
torque_profile(:,1) = fit_record(:,1);
%1st col: angle
%2nd col: average torque calculated from fit_record
%3rd col: average torque calculated from fit_record_lim
%4th col: average torque calculated from fit_record_rawdwshift
%5th col: average torque calculated from fit_record_rawdwshift_lim
%6th col: average of 2nd - 5th col
torque_pos = NaN*ones(length(fit_record),7);
%positive torque
%convention like torque_profile, but for positive torque
torque_pos(:,1) = fit_record(:,1);
torque_neg = NaN*ones(length(fit_record),7);
%negative torque
%convention like torque_profile, but for negative torque
torque_neg(:,1) = fit_record(:,1);

%store standard deviation of torque
torque_std = NaN*ones(length(fit_record),13);
torque_std(:,1) = fit_record(:,1);

for i = 1 : length(fit_record)
    
    
    torque_profile(i,2) = friction*fit_record(i,8)*1000*pi/180;
    torque_profile(i,3) = friction*fit_record_lim(i,8)*1000*pi/180;
    torque_profile(i,4) = friction*fit_record_rawdwshift(i,8)*1000*pi/180;
    torque_profile(i,5) = friction*fit_record_rawdwshift_lim(i,8)*1000*pi/180;
    torque_profile(i,6) = mean(torque_profile(i,2:5));

    torque_neg(i,2) = friction*fit_record(i,2)*1000*pi/180;
    torque_neg(i,3) = friction*fit_record_lim(i,2)*1000*pi/180;
    torque_neg(i,4) = friction*fit_record_rawdwshift(i,2)*1000*pi/180;
    torque_neg(i,5) = friction*fit_record_rawdwshift_lim(i,2)*1000*pi/180;
    torque_neg(i,6) = mean(torque_neg(i,2:5));
    
    torque_pos(i,2) = friction*fit_record(i,3)*1000*pi/180;
    torque_pos(i,3) = friction*fit_record_lim(i,3)*1000*pi/180;
    torque_pos(i,4) = friction*fit_record_rawdwshift(i,3)*1000*pi/180;
    torque_pos(i,5) = friction*fit_record_rawdwshift_lim(i,3)*1000*pi/180;
    torque_pos(i,6) = mean(torque_pos(i,2:5));
    
    torque_std(i,2) = friction*fit_record(i,4)*pi*1000/180; %error of neg torque 
    torque_std(i,3) = friction*fit_record(i,5)*pi*1000/180; %error of pos torque 
    
    torque_std(i,4) = friction*fit_record_lim(i,4)*pi*1000/180;
    torque_std(i,5) = friction*fit_record_lim(i,5)*pi*1000/180;
    
    torque_std(i,6) = friction*fit_record_rawdwshift(i,4)*pi*1000/180;
    torque_std(i,7) = friction*fit_record_rawdwshift(i,5)*pi*1000/180;
    
    torque_std(i,8) = friction*fit_record_rawdwshift_lim(i,4)*pi*1000/180;
    torque_std(i,9) = friction*fit_record_rawdwshift_lim(i,5)*pi*1000/180;
    
    
    if fit_record(i,6) == 1
        torque_std(i,10) = torque_std(i,2); 
                
        torque_std(i,11) = torque_std(i,4); 
                
        torque_std(i,12) = torque_std(i,6); 
               
        torque_std(i,13) = torque_std(i,8); 
                
    else        
        torque_std(i,10) = sqrt(sum((fit_record(i,6:7).*torque_std(i,2:3)).^2));
        %sqrt((fit_record(i,4)*torque_std(i,2))^2+(fit_record(i,5)*torque_std(i,3))^2); 
        torque_std(i,11) = sqrt(sum((fit_record_lim(i,6:7).*torque_std(i,2:3)).^2));
        torque_std(i,12) = sqrt(sum((fit_record_rawdwshift(i,6:7).*torque_std(i,2:3)).^2));
        torque_std(i,13) = sqrt(sum((fit_record_rawdwshift_lim(i,6:7).*torque_std(i,2:3)).^2));
    end
    
    torque_std(i,14) = sqrt(sum(torque_std(i,2:2:8).^2))/4;
    torque_std(i,15) = sqrt(sum(torque_std(i,3:2:9).^2))/4;
    torque_std(i,16) = sqrt(sum(torque_std(i,10:13).^2))/4;
end

torque_smooth = smooth(torque_profile(:,1),torque_profile(:,6),5,'moving');
torquepos_smooth= smooth(torque_pos(:,1),torque_pos(:,6),5,'moving');
torqueneg_smooth= smooth(torque_neg(:,1),torque_neg(:,6),5,'moving');

%plot torque for each data set 
figure
plot(fit_record(:,1), friction*fit_record(:,8)*1000*pi/180,'k.-')
hold on
plot(fit_record(temp,1), friction*fit_record(temp,2)*1000*pi/180,'b.-')
plot(fit_record(temp,1), friction*fit_record(temp,3)*1000*pi/180,'r.-')

plot(fit_record_lim(:,1), friction*fit_record_lim(:,8)*1000*pi/180,'k.-')
hold on
plot(fit_record_lim(temp,1), friction*fit_record_lim(temp,2)*1000*pi/180,'b.-')
plot(fit_record_lim(temp,1), friction*fit_record_lim(temp,3)*1000*pi/180,'r.-')

plot(fit_record_rawdwshift(:,1), friction*fit_record_rawdwshift(:,8)*1000*pi/180,'k.-')
hold on
plot(fit_record_rawdwshift(temp,1), friction*fit_record_rawdwshift(temp,2)*1000*pi/180,'b.-')
plot(fit_record_rawdwshift(temp,1), friction*fit_record_rawdwshift(temp,3)*1000*pi/180,'r.-')

plot(fit_record_rawdwshift_lim(:,1), friction*fit_record_rawdwshift_lim(:,8)*1000*pi/180,'k.-')
hold on
plot(fit_record_rawdwshift_lim(temp,1), friction*fit_record_rawdwshift_lim(temp,2)*1000*pi/180,'b.-')
plot(fit_record_rawdwshift_lim(temp,1), friction*fit_record_rawdwshift_lim(temp,3)*1000*pi/180,'r.-')

xlim([-4 121])
xlabel('Angle (deg)')
ylabel('Torque (pN.nm)')

%plot torque of each data set with error bar 
figure
errorbar(fit_record(:,1), friction*fit_record(:,8)*1000*pi/180,torque_std(:,10),'k.-')
hold on
errorbar(fit_record(temp,1), friction*fit_record(temp,2)*1000*pi/180,torque_std(temp,2),'b.-')
errorbar(fit_record(temp,1), friction*fit_record(temp,3)*1000*pi/180,torque_std(temp,3),'r.-')

errorbar(fit_record_lim(:,1), friction*fit_record_lim(:,8)*1000*pi/180,torque_std(:,11),'k.-')
hold on
errorbar(fit_record_lim(temp,1), friction*fit_record_lim(temp,2)*1000*pi/180,torque_std(temp,4),'b.-')
errorbar(fit_record_lim(temp,1), friction*fit_record_lim(temp,3)*1000*pi/180,torque_std(temp,5),'r.-')

errorbar(fit_record_rawdwshift(:,1), friction*fit_record_rawdwshift(:,8)*1000*pi/180,torque_std(:,12),'k.-')
hold on
errorbar(fit_record_rawdwshift(temp,1), friction*fit_record_rawdwshift(temp,2)*1000*pi/180,torque_std(temp,6),'b.-')
errorbar(fit_record_rawdwshift(temp,1), friction*fit_record_rawdwshift(temp,3)*1000*pi/180,torque_std(temp,7),'r.-')

errorbar(fit_record_rawdwshift_lim(:,1), friction*fit_record_rawdwshift_lim(:,8)*1000*pi/180,torque_std(:,13),'k.-')
hold on
errorbar(fit_record_rawdwshift_lim(temp,1), friction*fit_record_rawdwshift_lim(temp,2)*1000*pi/180,torque_std(temp,8),'b.-')
errorbar(fit_record_rawdwshift_lim(temp,1), friction*fit_record_rawdwshift_lim(temp,3)*1000*pi/180,torque_std(temp,9),'r.-')

xlim([-4 121])
xlabel('Angle (deg)')
ylabel('Torque (pN.nm)')

%plot smooth torque profile with error bar
figure
errorbar(torque_profile(:,1), torque_smooth, torque_std(:,16),'k.-')
hold on
errorbar(torque_pos(temp,1), torquepos_smooth(temp), torque_std(temp,15),'r.-')
errorbar(torque_neg(temp,1), torqueneg_smooth(temp), torque_std(temp,14),'b.-')
xlim([-4 121])
xlabel('Angle (deg)')
ylabel('Torque (pN.nm)')
legend('Average','Positive','Negative')

%% torque analysis by 2 methods
%method 1: linear fitting based on data of histogram
%method 2: from Gaussian fitting to have torque = 2mu/(kBT*sigma^2)
%use data obtained from fitting Gaussian of angvel rather than rerun
%--> for reproducibility, avoid random number generator by matlab in GMM 

torque = NaN*ones(length(fit_record),5);
%1st col: angular position
%2nd col: negative torque
%3rd col: positive torque
%4th col: average torque 
%for 1 gaussian, 2nd col = 3rd col = 4th col
%5th col: from linear fit of data based on fluctuation theorem 

%method 2: torque obtained from Gaussian fitting 
for i = 1 : length(torque)
    torque(i,1) = fit_record(i,1);
    %Convert angular velocity and standard deviation from deg/ms to rad 
    torque(i,2) = 2*4.14*(fit_record(i,2)*pi/(180*100))/(fit_record(i,4)*pi/(180*100))^2;   %1kT = 4.14 pN.nm
    torque(i,3) = 2*4.14*(fit_record(i,3)*pi/(180*100))/(fit_record(i,5)*pi/(180*100))^2;    %1kT = 4.14 pN.nm
    if fit_record(i,6) == fit_record(i,7)
        torque(i,4) = torque(i,2).*fit_record(i,6); %+ torque(i,3)*fit_record(i,7)    
    else
        torque(i,4) = sum(torque(i,2:3).*fit_record(i,6:7)); %+ torque(i,3)*fit_record(i,7)    
    end
end

%plot torque based on method 2
plot(torque(:,1),torque(:,2),'b.-',torque(:,1),torque(:,3),'r.-')
hold on
plot(torque(:,1),torque(:,4),'k.-')
legend('Negative torque','Positive torque','Average torque')
xlabel('Angular position (deg)')
ylabel('Torque (pN.nm)');
title('Average torque - Gaussian fitting')

%method 1 

%check the torque by plotting with data based on fluctuation theorem 
%for laptop
%fpath = 'D:\Onedrive-NTU\OneDrive - Nanyang Technological University\PhD\rev analyze\Results\Torque';
%for desktop at school 
fpath = 'D:\Onedrive-Home\OneDrive - Nanyang Technological University\PhD\rev analyze\Results\Torque'; 

%for i = 1 : length(torque_data)
%    torqueFT = torque_data{i,2};
for i = 1 : length(hist_record)
    avel = hist_record{i,2};
    cvel = hist_record{i,3};
    torqueFT = torque_check(cvel, avel);
    %linear fit from linear regression 
    lm2 = fitlm(torqueFT(:,1), torqueFT(:,2),'linear');
    coeff = lm2.Coefficients.Estimate;  
    %1st row: intercept
    %2nd row: slope
    %curve = fit(torqueFT(:,1), torqueFT(:,2), 'poly1');
    %coeff2=coeffvalues(curve);
    %torque(i,5) = coeff2(1);
    torque(i,5) = coeff(2);
    
    g = figure
    scatter(torqueFT(:,1),torqueFT(:,2));
    hold on     
    plot(torqueFT(:,1), coeff(2)*torqueFT(:,1) + coeff(1),'k');      
    plot(torqueFT(:,1), torque(i,4)*torqueFT(:,1),'r')        
    legend('Actual data','Linear fit from data','From Gaussian fit','Location', 'northwest')
    xlabel('\Delta\theta/kT (rad/pN.nm)')
    ylabel('ln(p(\Delta\theta)/p(-\Delta\theta)');
    title(sprintf('Torque extraction at %d deg', hist_record{i,1}));
    saveas(gca, fullfile(fpath, sprintf('torque_angle = %d',hist_record{i,1})), 'png'); 
    saveas(gca, fullfile(fpath, sprintf('torque_angle = %d',hist_record{i,1})), 'fig'); 
    %title(sprintf('Torque extraction at %d deg', torque_data{i,1}));
    %saveas(gca, fullfile(fpath, sprintf('torque_angle = %d',torque_data{i,1})), 'png'); 
    %aveas(gca, fullfile(fpath, sprintf('torque_angle = %d',torque_data{i,1})), 'fig'); 
end
close all

%compare between 2 fittings
plot(torque(:,1),torque(:,4),'r.-')
hold on
plot(torque(:,1),torque(:,5),'b.-')
legend('Gaussian fit','Linear fit')
xlabel('Angular position (deg)')
ylabel('Torque (pN.nm)');
title('Average torque - angular position')
xlim([-4 121])

%% torque analysis by Sandor 
fpath = 'D:\Onedrive-NTU\OneDrive - Nanyang Technological University\PhD\rev analyze\Results\Torque\';

torque = NaN*ones(length(fit_record),5);

for i = 1 : length(hist_record)
    avel = hist_record{i,2};
    cvel = hist_record{i,3};
    torqueFT = torque_check(cvel, avel);
    agl_pos = hist_record{i,1};
    torque(i,1) = agl_pos;

    %[Nb,Nf,pb,pf] = torque_SD(torqueFT, fpath, agl_pos);
    %torque(i,2:5) =  [Nb,Nf,pb,pf];      
    
    [Nb,Nf,p_ratio] = torque_SD2(torqueFT, fpath, agl_pos);
    torque(i,2:4) = [Nb,Nf,p_ratio];
end


%figure for 4 parameters 
plot(torque(15:27,1), torque(15:27,3),'.-')
hold on
plot(torque(15:27,1), torque(15:27,2),'.-')
plot(torque(15:27,1), torque(15:27,2).*torque(15:27,4) + torque(15:27,3).*torque(15:27,5),'.-')

figure
plot(torque(:,1), torque(:,3),'r.-')
hold on
plot(torque(:,1), torque(:,2),'b.-')
plot(torque(:,1), torque(:,2).*torque(:,4) + torque(:,3).*torque(:,5),'k.-')
xlabel('Angle (deg)')
ylabel('Torque (pN.nm)')
legend('Positive','Negative','Average')

figure
plot(torque(:,1), torque(:,4),'b.-')
hold on
plot(torque(:,1), torque(:,5),'r.-')
xlabel('Angle (deg)')
ylabel('Proability')
legend('Negative','Positive')

%figure for 3 parameters
figure
plot(torque(:,1), torque(:,3),'r.-')
hold on
plot(torque(:,1), torque(:,2),'b.-')
xlabel('Angle (deg)')
ylabel('Torque (pN.nm)')
legend('Positive','Negative')

figure
plot(torque(:,1), torque(:,4),'.-')
xlabel('Angle (deg)')
ylabel('Probability ratio')


%% check free energy usage by Pi release
% 
%also 

xdata = torque(:,1);
ydata = [torque(:,4), torque(:,5)];  %from Gaussian fit, linear fit 

index = find(xdata < 40);

%calculate energy as area under the curve by trapezoidal rule
energy = [0,0];
[nrow ncol] = size(ydata);

for j = 1 : ncol
    ydata_interest = ydata(:,j);
    
    for i = 1 : length(index) - 1
        temp = 0.5*(ydata_interest(index(i))+ydata_interest(index(i+1)))*3*pi/180;
        energy(1,j) = energy(1,j) + temp;
    end
end

%% Obtain statistic of 2 subsequent steps from each current angles
agl_pos = -3;
step = 3;

index_store = cell(2,2);
%1st col: angle
%2nd - 5th col: store count_pattern: row index where the stat occurs for up-up, up-down, down-up and down-down

stat_store = NaN*ones(3,9);
%1st col: angle
%2nd - 5th col: probability for up-up, up-down, down-up and down-down
%6th - 9th col: conditional propbability for up , down given up; up, down for given down respectively
i = 1;
while agl_pos <= 130
    [count_pattern, number, cond_prop] = updownstat(d_ang_raw, agl_pos);
    %[count_pattern, number, cond_prop] = updownstat(d_ang_raw_dwell_shift, agl_pos);
    index_store{i,1} = agl_pos;
    index_store(i,2:5) = count_pattern;
    stat_store(i,1) = agl_pos;
    stat_store(i,2:5) = number;
    stat_store(i,6:9) = cond_prop;
    
    i = i+1;
    agl_pos = agl_pos + step; 
end
close all

for k = 2:5
    plot(stat_store(:,1), stat_store(:,k),'.-')
    hold on
end
xlabel('Angular positions')
ylabel('Percentage')
title('Probability for next 2 steps')
legend('up-up','up-down','down-up','down-down')

name = {'up-up','up-down','down-up','down-down'};
figure
for k = 2:5
    subplot(2,2,k-1)
    plot(stat_store(:,1), stat_store(:,k),'.-')
    xlabel('Angular positions')
    ylabel('Percentage')
    title(name{k-1})
    xlim([-5 135])
end

%Plot conditional probability 


figure
for k = 2:5
    %subplot(2,2,k-1)
    g = figure
    plot(stat_store(:,1), stat_store(:,k+4),'.-')
    xlabel('Angular positions')
    ylabel('Percentage')
    title(name{k-1})
    xlim([-5 135])
end

%% Plot angular velocity vs time
jump_dur = NaN*ones(length(selected),1);
for i = 1 : length(selected)
    jump_dur(i,1) =  length(find(d_ang_raw(:,3) == selected(i,1)));
end
hist(jump_dur); %statistic of transition time 
ax = gca;
ax.XTickLabel =ax.XTick*0.01;
xlabel('Duration (ms)')
ylabel('Counts')
title('Transition time between dwells');

dur = 14;
idx = find(jump_dur == dur);
time = (0:0.01:0.01*(dur-1))';
%plot all angular position vs time 
for i = 1 : length(idx)-1
    temp = find(d_ang_raw(:,3) == selected(idx(i,1),1));
    %plot(d_ang_raw(temp,2), d_ang_raw(temp,1), '.-')
    %figure
    subplot(5,5,i)
    plot(time, d_ang_raw(temp,2), '.-')    
end
%plot all angular velocity vs time 
for i = 1 : length(idx)-1
    temp = find(d_ang_raw(:,3) == selected(idx(i,1),1));
    %plot(d_ang_raw(temp,2), d_ang_raw(temp,1), '.-')
    %figure
    subplot(5,5,i)
    plot(time, d_ang_raw(temp,1), '.-')    
    ax = gca;
    ax.XTick = 0:0.01:dur*0.01;
    ax.Box = 'off';
    xlabel('time (ms)');
    ylabel('Angular velocity (deg/ms)');
    %ntick = length(ax.XTick);
    ax2 = axes('Position',ax.Position,'XAxisLocation','top', 'YAxisLocation','left','Color','none','YTick',[],'YTickLabel',[])
    %line(d_ang_raw(temp,2), d_ang_raw(temp,1), 'Parent',ax2)
    data2_tick = round(d_ang_raw(temp,2)); 
    ax2.XLim = ax.XLim;
    ax2.XTick = ax.XTick;
    ax2.XTickLabel = data2_tick;
    %data2_tick=linspace(round(d_ang_raw(temp(1),2)),round(d_ang_raw(temp(end),2)),ntick);
    %set(ax1,'ytick',round(data1_tick*100)/100)
    %ax2.XLim = [round(d_ang_raw(temp(1),2)) round(d_ang_raw(temp(end),2))];
    %set(ax2,'xtick',round(data2_tick))
    xlabel('Angle (deg)');
    title(sprintf('i = %d',i));
end

posvel = cell(4,5);
%1st col: jumping duration
%2nd col: summary of all relevant trajectories with similar jumping duration
%3rd col: angular position 
%4th col: median of angular velocity
%5th col: median of angular position 

for dur = 10:13
    idx = find(jump_dur == dur);
    time = (0:0.01:0.01*(dur-1))';
    traj = [];
    agl = [];
    for i = 1 : length(idx)
        temp = find(d_ang_raw(:,3) == selected(idx(i,1),1));
        traj = [traj,d_ang_raw(temp,1)];      
        agl = [agl, d_ang_raw(temp,2)];  
        
    end
    med_traj = NaN*ones(length(time),1);
    med_agl = NaN*ones(length(time),1);
    for i = 1 : length(time)
        med_traj(i,1) = median(traj(i,:));
        med_agl(i,1) = median(agl(i,:));
    end
    
    posvel{dur-9,1} = dur;
    posvel{dur-9,2} = traj;
    posvel{dur-9,3} = agl;
    posvel{dur-9,4} = med_traj;
    posvel{dur-9,5} = med_agl;
end

%plot boxplot of angular velocity distribution 
figure
for i = 1 : 4
   subplot(2,2,i) 
   boxplot(posvel{i,2}');  
   ax = gca;
   ax.XTickLabel = 0:0.01:posvel{i,1}*0.01;
   xlabel('Time (ms)');
   ylabel('Angular velocity (deg/ms)')
   %ylabel('Angle (deg)');
end

%plot boxplot of angular position distribution 
figure
for i = 1 : 4
   subplot(2,2,i) 
   boxplot(posvel{i,3}');  
   ax = gca;
   ax.XTickLabel = 0:0.01:posvel{i,1}*0.01;
   xlabel('Time (ms)');   
   ylabel('Angle (deg)');
end


%plot median of angular position 
figure
for i = 1 : 4
    plot(1:length(posvel{i,5}),posvel{i,5},'.-')
    hold on
end
xlabel('Time (ms)');
ylabel('Angle (deg)');
title('Median angle vs time')

%angular velocity selected 
%dur 10
idx_choose = [2,6,9,15];

%dur 11, 
idx_choose =  [4,8,10,11,12,13,14,16,17,18];

%dur 12
idx_choose = [1,4,5,6,9,10,12,13,15,16,19,21,22,25];

%dur 13 
idx_choose =  [2,10,11,14,17,18,19,20,22,23];

for i = 1 : length(idx_choose)
    temp = find(d_ang_raw(:,3) == selected(idx(idx_choose(i),1),1));
    %plot(d_ang_raw(temp,2), d_ang_raw(temp,1), '.-')
    figure
    plot(time, d_ang_raw(temp,1), '.-')    
    ax = gca;
    ax.XTick = 0:0.01:dur*0.01;
    ax.Box = 'off';
    xlabel('time (ms)');
    ylabel('Angular velocity (deg/ms)');
    %ntick = length(ax.XTick);
    ax2 = axes('Position',ax.Position,'XAxisLocation','top', 'YAxisLocation','left','Color','none','YTick',[],'YTickLabel',[])
    %line(d_ang_raw(temp,2), d_ang_raw(temp,1), 'Parent',ax2)
    data2_tick = round(d_ang_raw(temp,2)); 
    ax2.XLim = ax.XLim;
    ax2.XTick = ax.XTick;
    ax2.XTickLabel = data2_tick;
    %data2_tick=linspace(round(d_ang_raw(temp(1),2)),round(d_ang_raw(temp(end),2)),ntick);
    %set(ax1,'ytick',round(data1_tick*100)/100)
    %ax2.XLim = [round(d_ang_raw(temp(1),2)) round(d_ang_raw(temp(end),2))];
    %set(ax2,'xtick',round(data2_tick))
    xlabel('Angle (deg)');
    title(sprintf('i = %d',i));
    
end

traj = [];
agl = [];
for i = 1 : length(idx_choose)
    temp = find(d_ang_raw(:,3) == selected(idx(idx_choose(i),1),1));
    traj = [traj,d_ang_raw(temp,1)];      
    agl = [agl, d_ang_raw(temp,2)];
    plot(time, d_ang_raw(temp,2),'.-')
    ax = gca;
    ax.XTick = 0:0.01:12*0.01;
    ax.Box = 'off';
    xlabel('time (ms)');
    ylabel('Angular velocity (deg/ms)');
    hold on
    %plot(time, d_ang_raw(temp,1),'.-')

end

%% pick trajectories for illustration 
%idx 10
idx_choose = [2,6,15];
%idx = 11
idx_choose = [17,18];
%idx = 12
idx_choose = [16,19];
%idx = 13
idx_choose = [1,11];
%idx = 14
idx_choose = [3,6,9,19,28];
%idx 15
idx_choose = [10,11,15];


dur = 14; 
idx = find(jump_dur == dur);
time = (0:0.01:0.01*(dur-1))';

for i = 1 : length(idx_choose)
    temp = find(d_ang_raw(:,3) == selected(idx(idx_choose(i),1),1));
    plot(time, d_ang_raw(temp,2),'.-')
    hold on
end


%the most probable angle for next jump
agl_most = NaN*ones(length(fit_record),4);
%1st col: current angle
%for 1 torque:
%2nd - 4th col: most probable next angle
%for 2 torques: 
%2nd - 4th col: most probable next angle of going back, forward, and average 

agl_most(:,1) = fit_record(:,1);
agl_most(:,2) =  fit_record(:,1) + fit_record(:,2)/100;
agl_most(:,3) =  fit_record(:,1) + fit_record(:,3)/100;
agl_most(:,4) =  fit_record(:,1) + fit_record(:,8)/100;

figure
plot(agl_most(:,1), agl_most(:,4), '.-')
xlim([-5 130]);
xlabel('Angle (deg)')
ylabel('Angle (deg)')
title('Next most probable angle')

%obtain a hypothetical most probable path 
start_agl = 0; 
traj_most = NaN*ones(3,2);
traj_most(1,1) = 0;
traj_most(1,2) = start_agl;
i = 1;
while start_agl <= 114
   [val, idx] = min(abs(start_agl - agl_most(:,1)));    

   i = i + 1;
   traj_most(i,2) = agl_most(idx,4); 
   traj_most(i,1) = 0.01*(i-1);
   start_agl = traj_most(i,2); 
    
end

figure
plot(traj_most(:,1), traj_most(:,2), 'r.-', 'LineWidth', 2)
legend('Average trajectory','Location', 'SouthEast')
hold on
%pick item from %pick trajectories for illustration 

picklist = NaN*ones(5,2);
picklist(:,1) = [10;12;13;14;15];
picklist(:,2) = [2;1;19;3;11];
for i = 1 : length(picklist)
    dur = picklist(i,1); 
    idx = find(jump_dur == dur);
    time = (0:0.01:0.01*(dur-1))';
    temp = find(d_ang_raw(:,3) == selected(idx(picklist(i,2),1),1));
    plot(time, d_ang_raw(temp,2),'.-')
    hold on
end
xlabel('Time (ms)');
ylabel('Angle (deg)');

%% other time relation stuff
selected = unique(d_ang_rawextended(:,3));

selected = unique(d_ang_raw(:,3));

nstep = 35;
v = NaN*ones(length(selected),nstep);
agl = NaN*ones(length(selected),nstep);
for j = 1 : nstep
%     v_temp = [];
%     agl_temp = [];
    v_temp = NaN*ones(length(selected),1);
    agl_temp = NaN*ones(length(selected),1);

    for i = 1 : length(selected)
%        temp = find(d_ang_rawextended(:,3) == selected(i,1));     
         temp = find(d_ang_raw(:,3) == selected(i,1));     

       if length(temp) >= j
%         v_temp = [v_temp,d_ang_raw(temp(j),1)];        
%         agl_temp = [agl_temp,d_ang_raw(temp(j),2)];    
%           v_temp(i) = d_ang_rawextended(temp(j),1);
%           agl_temp(i) = d_ang_rawextended(temp(j),2);
          v_temp(i) = d_ang_raw(temp(j),1);
          agl_temp(i) = d_ang_raw(temp(j),2);
       else 
           continue
       end
   end
    v(:,j) = v_temp;
    agl(:,j) = agl_temp;
end

vmean = [];
aglmean = [];
for i = 1 : nstep
    vmean = [vmean, mean(v{i})];
    aglmean = [aglmean, mean(agl{i})];
    %plot(i,median(v{i}),'-kx')
    %hold on
end
nbin = 30;
count_list = NaN*ones(nstep,nbin);
edge = -10:5:130;
for i = 10: nstep
%    [count, val] = hist(agl(:,i),20); 
%    count_list(:,i) = count';
%      [count, val] = 
     hist(agl(:,i),nbin);
%      [count, val] = histogram(agl(:,i),edge);

%      count_list(i,:) = count;
     
     title(sprintf('step %d',i))
     xlabel('Angle (deg)')
     ylabel('Count')
     pause
end


bar3(count_list)

bar3(1:nstep, count_list')

plot(1:length(vmean),vmean,'-x')
plot(1:length(vmean),aglmean,'-x')
hist(agl{6},30)

boxplot(1:2,[v{1},v{2}])
 
 