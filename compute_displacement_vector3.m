function U = compute_displacement_vector3(Ja, Raj, d_a, alph, k, NN)
% COMPUTE_DISPLACEMENT_VECTOR  Deterministic SH-wave scattering solver.
%   U = compute_displacement_vector(xx, Ja, Raj, d_a, alph, k, NN)
%   computes the amplitude of surface displacement |U| for a group of
%   semi-cylindrical canyons under SH-wave incidence using wave function
%   expansion and Graf's addition theorem.
%
%   Inputs:
%       xx   - vector of observation point coordinates (normalized by a0)
%       Ja   - vector of canyon indices (e.g., [-2 -1 0 1 2])
%       Raj  - vector of canyon radii (normalized by a0)
%       d_a  - vector of distances from the central canyon (normalized by a0)
%       alph - incident angle (radians)
%       k    - wave number
%       NN   - truncation order for wave function expansion
%
%   Output:
%       U    - displacement amplitude at surface points (column vector)

    a0 = 1;                         % 基准半宽（与主程序一致）
    d = d_a * a0;                   % 实际距离
    NJ = length(Ja);
    xx = (-8:0.02:8)*a0;
    Numx = length(xx);
    
    % --- 构建并求解线性系统 ---
    Ln = zeros(NJ*NN, NJ*NN);
    Rn = zeros(NJ*NN, 1);
    
    for m = 0:NJ*NN-1
        Nump = floor(m/NN) + 1;
        p = Ja(Nump);
        if p < 0
            Dp0 = -d_a(Nump);
        elseif p == 0
            Dp0 = 0;
        else
            Dp0 = d_a(Nump);
        end
        
        Jm = m - (Nump-1)*NN;
        Rn(m+1,1) = -exp(1i*k*Dp0*cos(alph)) * 4 * 1i^Jm * cos(Jm*alph);
        
        for n = 0:NJ*NN-1
            Numj = floor(n/NN) + 1;
            Jn = n - (Numj-1)*NN;
            if Jn == 0
                epsn = 1;
            else
                epsn = 2;
            end
            
            if Numj < Nump
                if Ja(Numj)*Ja(Nump) < 0
                    DJP = d(Numj) + d(Nump);
                else
                    DJP = abs(d(Numj) - d(Nump));
                end
                Kmn = besselh(Jn+Jm, 1, k*DJP) + (-1)^Jm * besselh(Jn-Jm, 1, k*DJP);
                Ln(m+1,n+1) = (-1)^Jm * Kmn;
                
            elseif Numj == Nump
                rp = Raj(Nump) * a0;
                Dbesh = Jn/rp * besselh(Jn, 1, k*rp) - k * besselh(Jn+1, 1, k*rp);
                Dbesj = Jn/rp * besselj(Jn,     k*rp) - k * besselj(Jn+1,     k*rp);
                if Jn == Jm
                    Ln(m+1,n+1) = (2 * Dbesh) / (Dbesj * epsn);
                else
                    Ln(m+1,n+1) = 0;
                end
                
            else % Numj > Nump
                if Ja(Numj)*Ja(Nump) < 0
                    DPJ = d(Numj) + d(Nump);
                else
                    DPJ = abs(d(Nump) - d(Numj));
                end
                Kmn = besselh(Jn+Jm, 1, k*DPJ) + (-1)^Jm * besselh(Jn-Jm, 1, k*DPJ);
                Ln(m+1,n+1) = (-1)^Jn * Kmn;
            end
        end
    end
    
    An = Ln \ Rn;
    An1 = An(1:NN);
    An2 = An(NN+1:2*NN);
    An3 = An(2*NN+1:3*NN);
    
    % --- 计算地表位移（向量化处理） ---
    wf  = zeros(Numx, 1);
    ws1 = zeros(Numx, 1);
    ws2 = zeros(Numx, 1);
    ws3 = zeros(Numx, 1);
    
    for m = 1:Numx
        x0 = xx(m);
        
        % 根据x0位置确定坐标计算（原代码逻辑，此处保持完全一致）
        if x0 <= -d(1)-Raj(1)
            y0 = 0; r0 = sqrt(x0^2+y0^2); a0_ang = acos(x0/r0);
            yN1 = 0; xN1 = x0+d(1); rN1 = sqrt(xN1^2+yN1^2); aN1 = acos(xN1/rN1);
            yP1 = 0; xP1 = x0-d(3); rP1 = sqrt(xP1^2+yP1^2); aP1 = acos(xP1/rP1);
            wf(m) = exp(1i*k*(x0*cos(alph)-y0*sin(alph))) + exp(1i*k*(x0*cos(alph)+y0*sin(alph)));
            sums1 = 0; sums2 = 0; sums3 = 0;
            for n = 0:NN-1
                sums1 = sums1 + An1(n+1)*besselh(n,1,k*rN1)*cos(n*aN1);
                sums2 = sums2 + An2(n+1)*besselh(n,1,k*r0)*cos(n*a0_ang);
                sums3 = sums3 + An3(n+1)*besselh(n,1,k*rP1)*cos(n*aP1);
            end
            ws1(m) = sums1; ws2(m) = sums2; ws3(m) = sums3;
            
        elseif x0 <= -d(1)+Raj(1)
            xN1 = x0+d(1); yN1 = sqrt(Raj(1)^2-xN1^2); rN1 = sqrt(xN1^2+yN1^2); aN1 = acos(xN1/rN1);
            y0  = yN1;   r0 = sqrt(x0^2+y0^2); a0_ang  = acos(x0/r0);
            yP1 = yN1;  xP1 = x0-d(3); rP1 = sqrt(xP1^2+yP1^2); aP1 = acos(xP1/rP1);
            wf(m) = exp(1i*k*(x0*cos(alph)-y0*sin(alph))) + exp(1i*k*(x0*cos(alph)+y0*sin(alph)));
            sums1 = 0; sums2 = 0; sums3 = 0;
            for n = 0:NN-1
                sums1 = sums1 + An1(n+1)*besselh(n,1,k*rN1)*cos(n*aN1);
                sums2 = sums2 + An2(n+1)*besselh(n,1,k*r0)*cos(n*a0_ang);
                sums3 = sums3 + An3(n+1)*besselh(n,1,k*rP1)*cos(n*aP1);
            end
            ws1(m) = sums1; ws2(m) = sums2; ws3(m) = sums3;
            
        elseif x0 <= d(2)-Raj(2)
            y0 = 0; r0 = sqrt(x0^2+y0^2); a0_ang = acos(x0/r0);
            yN1 = 0; xN1 = x0+d(1); rN1 = sqrt(xN1^2+yN1^2); aN1 = acos(xN1/rN1);
            yP1 = 0; xP1 = x0-d(3); rP1 = sqrt(xP1^2+yP1^2); aP1 = acos(xP1/rP1);
            wf(m) = exp(1i*k*(x0*cos(alph)-y0*sin(alph))) + exp(1i*k*(x0*cos(alph)+y0*sin(alph)));
            sums1 = 0; sums2 = 0; sums3 = 0;
            for n = 0:NN-1
                sums1 = sums1 + An1(n+1)*besselh(n,1,k*rN1)*cos(n*aN1);
                sums2 = sums2 + An2(n+1)*besselh(n,1,k*r0)*cos(n*a0_ang);
                sums3 = sums3 + An3(n+1)*besselh(n,1,k*rP1)*cos(n*aP1);
            end
            ws1(m) = sums1; ws2(m) = sums2; ws3(m) = sums3;
            
        elseif x0 <= d(2)+Raj(2)
            y0  = sqrt(Raj(2)^2-x0^2);  r0  = sqrt(x0^2+y0^2); a0_ang  = acos(x0/r0);
            xN1 = x0+d(1); yN1 = y0; rN1 = sqrt(xN1^2+yN1^2); aN1 = acos(xN1/rN1);
            xP1 = x0-d(3); yP1 = y0; rP1 = sqrt(xP1^2+yP1^2); aP1 = acos(xP1/rP1);
            wf(m) = exp(1i*k*(x0*cos(alph)-y0*sin(alph))) + exp(1i*k*(x0*cos(alph)+y0*sin(alph)));
            sums1 = 0; sums2 = 0; sums3 = 0;
            for n = 0:NN-1
                sums1 = sums1 + An1(n+1)*besselh(n,1,k*rN1)*cos(n*aN1);
                sums2 = sums2 + An2(n+1)*besselh(n,1,k*r0)*cos(n*a0_ang);
                sums3 = sums3 + An3(n+1)*besselh(n,1,k*rP1)*cos(n*aP1);
            end
            ws1(m) = sums1; ws2(m) = sums2; ws3(m) = sums3;
            
        elseif x0 <= d(3)-Raj(3)
            y0  = 0; r0  = sqrt(x0^2+y0^2);   a0_ang = acos(x0/r0);
            yN1 = 0; xN1 = x0+d(1); rN1 = sqrt(xN1^2+yN1^2); aN1 = acos(xN1/rN1);
            yP1 = 0; xP1 = x0-d(3); rP1 = sqrt(xP1^2+yP1^2); aP1 = acos(xP1/rP1);
            wf(m) = exp(1i*k*(x0*cos(alph)-y0*sin(alph))) + exp(1i*k*(x0*cos(alph)+y0*sin(alph)));
            sums1 = 0; sums2 = 0; sums3 = 0;
            for n = 0:NN-1
                sums1 = sums1 + An1(n+1)*besselh(n,1,k*rN1)*cos(n*aN1);
                sums2 = sums2 + An2(n+1)*besselh(n,1,k*r0)*cos(n*a0_ang);
                sums3 = sums3 + An3(n+1)*besselh(n,1,k*rP1)*cos(n*aP1);
            end
            ws1(m) = sums1; ws2(m) = sums2; ws3(m) = sums3;
            
        elseif x0 <= d(3)+Raj(3)
            xP1 = x0-d(3); yP1 = sqrt(Raj(3)^2-xP1^2); rP1 = sqrt(xP1^2+yP1^2); aP1 = acos(xP1/rP1);
            y0 = yP1; r0 = sqrt(x0^2+y0^2); a0_ang = acos(x0/r0);
            xN1 = x0+d(1); yN1 = yP1; rN1 = sqrt(xN1^2+yN1^2); aN1 = acos(xN1/rN1);
            wf(m) = exp(1i*k*(x0*cos(alph)-y0*sin(alph))) + exp(1i*k*(x0*cos(alph)+y0*sin(alph)));
            sums1 = 0; sums2 = 0; sums3 = 0;
            for n = 0:NN-1
                sums1 = sums1 + An1(n+1)*besselh(n,1,k*rN1)*cos(n*aN1);
                sums2 = sums2 + An2(n+1)*besselh(n,1,k*r0)*cos(n*a0_ang);
                sums3 = sums3 + An3(n+1)*besselh(n,1,k*rP1)*cos(n*aP1);
            end
            ws1(m) = sums1; ws2(m) = sums2; ws3(m) = sums3;
            
        else
            y0 = 0; r0 = sqrt(x0^2+y0^2); a0_ang = acos(x0/r0);
            yN1 = 0; xN1 = x0+d(1); rN1 = sqrt(xN1^2+yN1^2); aN1 = acos(xN1/rN1);
            yP1 = 0; xP1 = x0-d(3); rP1 = sqrt(xP1^2+yP1^2); aP1 = acos(xP1/rP1);
            wf(m) = exp(1i*k*(x0*cos(alph)-y0*sin(alph))) + exp(1i*k*(x0*cos(alph)+y0*sin(alph)));
            sums1 = 0; sums2 = 0; sums3 = 0;
            for n = 0:NN-1
                sums1 = sums1 + An1(n+1)*besselh(n,1,k*rN1)*cos(n*aN1);
                sums2 = sums2 + An2(n+1)*besselh(n,1,k*r0)*cos(n*a0_ang);
                sums3 = sums3 + An3(n+1)*besselh(n,1,k*rP1)*cos(n*aP1);
            end
            ws1(m) = sums1; ws2(m) = sums2; ws3(m) = sums3;
        end
    end
    
    U = abs(wf + ws1 + ws2 + ws3);
end