function [min_error,ti] = inflectpt(time,perimeter)
    
    
    % Code meant to detect the point of inflexion of a curve provide
    % This code was developped by Alexis Villars in the Team of Romain Levayer at Institut Pasteur 
    % Last modification 09/05/2022

    % This code is provided under an MIT license in a repository for the publication in nature communications
    
    % This functions computes two fits. One between t0 and ti and one between ti and t_end. For each ti, it computes the error between the fit 
    % and the actual datas in perimeter. 
    
    % it return two variables
    % min_error is the minimum of the error found for all ti 
    % ti is the inflexion point obtained for the min_error. 
    
    % it takes 2 parameters
    % time is a vector of the same length of the perimeter vector. Best pratice is to set it to 1:length(perimeter) 
    % to return direct index of the found inflexion point
    % perimeter is the actual vector of data for which you want to found the inflexion point. 
    
    error_total = zeros(length(time)-3,1);
         
    for jj = 2:length(time)-1

        f1=fit(time(1:jj),perimeter(1:jj),'poly1');
        f2=fit(time(jj:end),perimeter(jj:end),'poly1');

        error_f1 = 0;
        for t = 1:jj
            error_1 = abs(f1(t)-perimeter(t));
            error_f1 = error_f1+error_1;       
        end        
        
        error_f2 = 0; 
        for t = jj:length(time)
            error_2 = abs(f2(t)-perimeter(t));
            error_f2 = error_f2+error_2;
        end
        
        error_final = error_f1+error_f2;
        error_total(jj-1,1) = error_final; 
        
    end
    
    [min_error,ti] = min(error_total);
    
end
