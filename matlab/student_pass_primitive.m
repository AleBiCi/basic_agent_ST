%#codegen
%%
%           Agent Logic
%         Pass Primitive
%              2025
%
%
%%

function [coeffsT2, v2, T2, coeffsT1, v1, T1] = student_pass_primitive(v0, a0, sf, vfmin, vfmax, Tmin, Tmax)
    if a0 >= 0
        Tvmin = final_opt_time_pass(v0, a0, sf, vfmin);
        Tvmax = final_opt_time_pass(v0, a0, sf, vfmax);
    else
        T_star = time_min_vel(a0, sf);
        v_star = min_vel(v0, a0, sf);
        % <= so that FreeFlow (vfmin == vfmax == vr, slide 104 line 3) still has a solution when a0 < 0
        if v_star < vfmin && vfmin <= vfmax
            Tvmin = final_opt_time_pass(v0, a0, sf, vfmin);
            Tvmax = final_opt_time_pass(v0, a0, sf, vfmax);
        elseif vfmin < v_star && v_star < vfmax
            Tvmin = T_star;
            Tvmax = final_opt_time_pass(v0, a0, sf, vfmax);
        else
            Tvmin = 0.;
            Tvmax = 0.;
        end
    end

    % FREE-FLOW PRIMITIVE
    if Tmin == 0. && Tmax == 0.
        T1 = Tvmax;     % fastest (reaches vfmax), slide 59 line 18 ordering
        T2 = Tvmin;
    % FREE-FLOW PRIMITIVE
    else
        T1 = max(Tmin, Tvmax);
        T2 = min(Tmax, Tvmin);
    end

    if Tvmax ~= 0 && Tvmax <= Tvmin && T1 > 0 && T1 <= T2
        v1 = final_opt_vel_pass(v0, a0, sf, T1);
        v2 = final_opt_vel_pass(v0, a0, sf, T2);
        coeffsT1 = coef_list_fun(v0, a0, sf, v1, 0., T1);
        coeffsT2 = coef_list_fun(v0, a0, sf, v2, 0., T2);
    else
        coeffsT1 = zeros(1,6);  % row, same size as coef_list_fun -> fixed-size double[6] in C
        coeffsT2 = zeros(1,6);
        T1 = 0.;
        T2 = 0.;
        v1 = 0.;
        v2 = 0.;
    end
end
