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

