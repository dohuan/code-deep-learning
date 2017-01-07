function [coeff, Ortho] = getCM(para_info, flag)
if flag==1 % --- load file
    M = load('./data/coeff');
    coeff = M.coeff;
    M = load('./data/Ortho');
    Ortho = M.Ortho;
else
    
    k_sigma_f = 0.05;
    k_sigma_m = 0.01;
    Le_lowH = 12.027458;             % patient H
    Le_lowK = 1.302589479770089e+01; % patient K
    
    run_day = 3000;
    np = size(para_info,2);
    poly_order = 3;
    for i=1:np
        [Ortho{i},para_value(:,i)] = coll_points_generate(poly_order,para_info(i));
    end
    coll_pts = [];
    for k1=1:poly_order+1
        for k2=1:poly_order+1
            for k3=1:poly_order+1
                coll_pts = [coll_pts;[para_value(k1,1) para_value(k2,2)...
                    para_value(k3,3)]];
            end
        end
    end
    y_fem = [];
    parfor i=1:size(coll_pts,1)
        % --- RHS
        %out = AAA_main(0.01, 0.01, run_day, 'test','false',coll_pts(i,:));
        %y_fem(i,1) = max(out.max_diameter);
        damage_k = coll_pts(i,1);
        damage_sigma = coll_pts(i,2);
        damage_mu = coll_pts(i,3);
        tmp = AAA_main(k_sigma_f,k_sigma_m,Le_lowH, damage_k,...
            damage_sigma,damage_mu, run_day,'H','true');
        y_fem = [y_fem; reshape(tmp',1,[])];
        
        % --- LHS
        H_temp = [];
        
        for k1=1:poly_order+1
            for k2=1:poly_order+1
                for k3=1:poly_order+1
                    H_temp = [H_temp,Ortho{1}.H{k1}(coll_pts(i,1))*Ortho{2}.H{k2}(coll_pts(i,2))*...
                            Ortho{3}.H{k3}(coll_pts(i,3))];
                end
            end
        end
        
        K(i,:) = H_temp;
        fprintf('Create CM equations... %d%%\n',round(i/size(coll_pts,1)*100));
    end
    %delete(poolobj);
    
    nancount = 0;
    mean_y = nanmean(y_fem,1);
    for i=1:size(y_fem,1)
        if (sum(isnan(y_fem(i,:)))~=0)
            y_fem(i,:) = mean_y;
            nancount = nancount + 1;
        end
    end
    fprintf('There are %d NaNs in y_fem\n',nancount);
    coeff = K\y_fem;
    
    save('./data/coeff','coeff');
    save('./data/Ortho','Ortho');
end
end