%#codegen
%%
%          Agent Logic
%         Pass Primitive J0
%              2025
%
%
%%

function [coefsj0, vfj0, tfj0] = student_pass_primitive_j0(v0, a0, sf, vfmin, vfmax)
    % Real roots only if 8 a0 sf + 5 v0^2 >= 0 (slide 90); a root is valid only if T > 0
    if v0 <= 0 || sf <= 0 || 8*a0*sf + 5*v0^2 < 0
        coefsj0 = zeros(1,6);
        tfj0 = 0.;
        vfj0 = 0.;
        return;
    end
    tfj0_list = final_opt_time_pass_j0(v0, a0, sf);
    vfj0 = final_opt_vel_pass_j0(v0,a0,0,sf,0,tfj0_list(1));
    if tfj0_list(1) > 0 && vfmin < vfj0 && vfmax > vfj0
        tfj0 = tfj0_list(1);
        coefsj0 = coef_list_fun(v0, a0, sf, vfj0, 0., tfj0);
    else
        vfj0 = final_opt_vel_pass_j0(v0,a0,0,sf,0,tfj0_list(2));
        if tfj0_list(2) > 0 && vfmin < vfj0 && vfmax > vfj0
            tfj0 = tfj0_list(2);
            coefsj0 = coef_list_fun(v0, a0, sf, vfj0, 0., tfj0);
        else
            coefsj0 = zeros(1,6);
            tfj0 = 0.;
            vfj0 = 0.;
        end
    end
end