%#codegen
%%
%           Agent Logic
%         Pass Primitive
%              2025
%
%
%%

function [coeffsT2, v2, T2, coeffsT1, v1, T1] = student_pass_primitive(v0, a0, sf, vfmin, vfmax, Tmin, Tmax)
    Tvmin = 0.;
    Tvmax = 0.;
    if a0 >= 0
        Tvmin = final_opt_time_pass(v0, a0, sf, vfmin);
        Tvmax = final_opt_time_pass(v0, a0, sf, vfmax);
    else
        
    end
end
