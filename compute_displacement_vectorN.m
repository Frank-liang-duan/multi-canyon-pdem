function U = compute_displacement_vectorN(Ja, Raj, d_a, alph, k, NN)
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

a0 = 1;                         
d = d_a * a0;                   
NJ = length(Ja);
xxxx = (-1.0*a0:0.02*a0:1.0*a0)';
Numx = length(xxxx);


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
An = Ln\Rn;

wf  = zeros(Numx,1);
ws  = wf;
CN  = (NJ+1)/2;

for m = 1:Numx
    x0 = xxxx(m);
    y0 = sqrt(Raj(CN)-x0^2);
    xc = d(CN)*Ja(CN);    
    xf = xc+x0; yf = y0;
    wf(m) = exp(1i*k*(xf*cos(alph)-yf*sin(alph)))+exp(1i*k*(xf*cos(alph)+yf*sin(alph)));
    sums = 0;
    for mj = 1:NJ
        Cn = An((mj-1)*NN+1:mj*NN);
        xj = d(mj)*Ja(mj);    
        if mj < CN           
            xx = (xc-xj)+x0; yy = y0;
            rN = sqrt(xx^2+yy^2); aN = acos(xx/rN);
            temp = 0;
            for n = 0:NN-1
                temp = temp+Cn(n+1)*besselh(n,1,k*rN)*cos(n*aN);
            end
        elseif mj == CN  
            xx = x0; yy = y0;
            rN = sqrt(xx^2+yy^2); aN = acos(xx/rN);
            temp = 0;
            for n = 0:NN-1
                temp = temp+Cn(n+1)*besselh(n,1,k*rN)*cos(n*aN);
            end
        elseif mj > CN   % 表示计算河谷右侧的河谷
            xx =-(xj-xc)+x0; yy = y0;
            rN = sqrt(xx^2+yy^2); aN = acos(xx/rN);
            temp = 0;
            for n = 0:NN-1
                temp = temp+Cn(n+1)*besselh(n,1,k*rN)*cos(n*aN);
            end
        end
        sums = sums+temp;
    end
    ws(m) = sums;

end

U   = abs(wf+ws);
end