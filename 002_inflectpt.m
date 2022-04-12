function [min_error,ti] = inflectpt(time,perimeter)

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