function [numComponents, MU, COV, PPp, histdata] = GMM(d_ang_raw120, agl_pos)
    %Use Gaussian mixture model to fit angular velocity distribution
    
    step = 3;
    epsilon = 3;
    bin = 20;

    temp = find( agl_pos -epsilon <= d_ang_raw120(:,2) & d_ang_raw120(:,2) <= agl_pos+epsilon);

    interest =  d_ang_raw120(temp(:,1),1);
    %for applying limit case 
    temp2 = find(interest < -5000 | interest > 6500);
    if isempty(temp2) ~= 1
        interest(temp2) = [];
    end
    
    [cvel, avel] = hist(interest,bin);
    
    %to eliminate extreme count due to uncertainty in knowing 
    %starting and ending of transition 
    if agl_pos < 10 | agl_pos > 110
        [value, index] = max(cvel);
        cvel(index) = [];
        avel(index) = [];
    end
   
    
    binwidth = avel(2) - avel(1);

    %find number of gaussian 
    AIC = zeros(1,2); BIC = zeros(1,2);
    obj = cell(1,2);
    options = statset('MaxIter',500); 
    for k = 1:2
       obj{k} = gmdistribution.fit( interest, k,'CovType','diagonal', 'Options',options);
       %obj{k} = gmdistribution.fit( d_ang(temp(:,1),1), k, 'Options',options);
       AIC(k) = obj{k}.AIC;    %based on AIC Akaike information
       BIC(k) = obj{k}.BIC;    %based on BIC Bayesian information
    end

    %[minAIC,numComponentsAIC] = min(AIC); 
    %[minBIC,numComponentsBIC] = min(BIC);
   
    histdata = {avel;cvel;obj{1};obj{2}}; 
    
    %fit based on number of gaussian 
    if   35 <= agl_pos && agl_pos <= 70
        numComponents = 2;
    else numComponents = 1;
    end
    
    %to fit at all angles 1 Gaussian 
    %numComponents = 1; 
    
    rng('default'); %to fix random generator --> for reproducibility 
    rng(1);
    paramEsts= gmdistribution.fit(interest,numComponents,'Replicates',1,'CovType','diagonal');  
    
    %paramEstsAIC = gmdistribution.fit(interest,numComponentsAIC);  
    %paramEstsBIC = gmdistribution.fit(interest,numComponentsBIC);  

    
    %histdata = {avel;cvel;minAIC;numComponentsAIC;paramEstsAIC;minBIC;numComponentsBIC;paramEstsBIC};   
    
    switch numComponents
        case 1
             MU = paramEsts.mu;
             COV = paramEsts.Sigma;   %sigma here is the covariance 
             PPp = paramEsts.PComponents;
             %fit_info(count,2) = MU;
             %fit_info(count,3) = SIGMA;
        case 2
            MU=[paramEsts.mu(1);paramEsts.mu(2);];
            COV = cat(3,[paramEsts.Sigma(1)],[paramEsts.Sigma(2)]);  %store covariance 
            PPp = [paramEsts.PComponents(1),paramEsts.PComponents(2)];  %component proportion
            %fit_info(count,2:3) = MU;
            %fit_info(count,4:5) = SIGMA;
            %fit_info(count,6:7) = PPp;
    end
    objA = gmdistribution(MU,COV,PPp);
    xgridss=transpose(linspace(min(avel),max(avel),100)); 
     
    %plot pdf 
    %for laptop
    fpath = 'D:\Onedrive-NTU\OneDrive - Nanyang Technological University\PhD\rev analyze\Results\Ang vel\';
    %for desktop at school 
    %fpath = 'D:\Onedrive-Home\OneDrive - Nanyang Technological University\PhD\rev analyze\Results\Ang vel'; 
    %fpath = [fpath,sprintf('%d',agl_pos)];
    g = figure
    %plot(avel,cvel/sum(cvel*binwidth),'ko'); hold on %convert to pdf so that area under curve = 1
    %plot(xgridss,pdf(objA,xgridss),'r-','linewidth',2)
    plot(avel,cvel,'ko'); 
    %histogram(interest,bin);
    hold on 
    plot(xgridss,pdf(objA,xgridss)*sum(cvel*binwidth),'r-','linewidth',2)  %convert pdf to real statistic    
    if numComponents == 2  
        %backward is blue, forward is orange
        if MU(1) < MU(2)
            obj_1 = gmdistribution(MU(1), COV(1), PPp(1));
            obj_2 = gmdistribution(MU(2), COV(2), PPp(2));
            plot(xgridss,pdf(obj_1,xgridss)*sum(cvel*binwidth)*PPp(1),'Color',[0 0.45 0.74])
            plot(xgridss,pdf(obj_2,xgridss)*sum(cvel*binwidth)*PPp(2),'Color',[0.85 0.33 0.1])
        else
            obj_2 = gmdistribution(MU(1), COV(1), PPp(1));
            obj_1 = gmdistribution(MU(2), COV(2), PPp(2));
            plot(xgridss,pdf(obj_1,xgridss)*sum(cvel*binwidth)*PPp(2),'Color',[0 0.45 0.74])
            plot(xgridss,pdf(obj_2,xgridss)*sum(cvel*binwidth)*PPp(1),'Color',[0.85 0.33 0.1])
        end
       
        %plot(xgridss,pdf(obj_1,xgridss)*sum(cvel*binwidth)*PPp(1) + pdf(obj_2,xgridss)*sum(cvel*binwidth)*PPp(2))
    end
    title(sprintf('Angular velocity distribution at %d deg',agl_pos))
    xlabel('Angular velocity (deg/ms)')
    ylabel('Counts')
    %saveas(gca, fullfile(fpath, sprintf('angle = %d',agl_pos)), 'png');   
    %saveas(gca, fullfile(fpath, sprintf('angle = %d',agl_pos)), 'fig'); 
    %close all
           
end
