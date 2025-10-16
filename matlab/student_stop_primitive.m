%#codegen
%%
%          Agent Logic
%         Stop Primitive
%              2025
%%

function [coefs,maxsf,tf] = student_stop_primitive(v0,a0,sf)
    % Impossible cases: negative velocity and null final position -->
    % return 0
    if (v0 <= 0) || (sf == 0)
        coefs = zeros(1,6); % Returns a vector of 6 zeros for each coefficient
        tf = 0.;
        maxsf = 0.; % Also set other returns to 0
    else
        % Can't reach a specified sf using given inputs (v0, a0)
        if 4*v0^2 + 5*a0*sf < 0
            maxsf = - (4*v0^2)/(5*a0);
            tf = (10*maxsf)/(2*v0);
        else
            % Valid case: compute final time with defined function
            maxsf = sf;
            tf = final_opt_time_stop(v0, a0, maxsf);
        end
        
        % Evaluate coefficients (consider vf = af = 0)
        coefs = coef_list_fun(v0, a0, 0., 0., 0., tf);
    end
end