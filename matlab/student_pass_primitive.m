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
    T_star = 0.;
    v_star = 0.;
    v1 = 0.;
    v2 = 0.;
    if a0 >= 0
        Tvmin = final_opt_time_pass(v0, a0, sf, vfmin);
        Tvmax = final_opt_time_pass(v0, a0, sf, vfmax);
    else
        T_star = time_min_vel(v0, a0, sf, T);
        v_star = min_vel(v0, a0, sf, T);
        if v_star < v_min && v_min < v_max
            Tvmin = final_opt_time_pass(v0, a0, sf, vfmin);
            Tvmax = final_opt_time_pass(v0, a0, sf, vfmax);
        elseif vfmin < v_star && v_star < vfmax
            Tvmin = T_star;
            Tvmax = final_opt_time_pass(v0, a0, sf, vfmax);
        end
    end
    T1 = 0.; T2 = 0.;
    [T1, T2] = intersect([Tmin, Tmax],[Tvmax, Tvmin]);

    if T1 > 0. & T1 <= T2
        v1 = final_opt_vel_pass(v0, a0, sf, T2);
        v2 = final_opt_vel_pass(v0, a0, sf, T1);
        coeffsT1 = coef_list_fun(v0, a0, sf, v1, 0., T1);
        coeffsT2 = coef_list_fun(v0, a0, sf, v2, 0., T2);
    end

    coeffsT1 = zeros(6);
    coeffsT2 = zeros(6);
end
