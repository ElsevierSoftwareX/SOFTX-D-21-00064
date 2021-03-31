function [r, ires] = res_dassl_dim (t, w, wd)
    
// --------------------------------------------------------- //
// RHS function                                              //
// Last modified: 02/04/2021                                 //
// --------------------------------------------------------- //
    
    r       = zeros(N3,1);  // Initialization
    dphi    = w(1:N);       // Electrode voltage [m^2.kg.s^-3.A^-1]
    eps     = w(N+1:N2);    // Porosity
    cA      = w(N2+1:N3);   // Acid concentration [mol.m^-3]
    dphit   = wd(1:N);      // Electrode voltage derived as a function of time [m^2.kg.s^-4.A^-1]
    epst    = wd(N+1:N2);   // Porosity derived as a function of time [s^-1]
    cAt     = wd(N2+1:N3);  // Acid concentration derived as a function of time [mol.s^-1.m^-3]


//------------------------ Physical relations ---------------------------------

    mA_val      = mA(cA);                       // Sulfuric acid molality [mol.kg^-1]
    muA_ref     = muA(mAref);                   // Reference electrochemical potential of sulfuric acid [m^2.kg.s^-2.mol^-1]
    muA_val     = muA(mA_val);                  // Electrochemical potential of sulfuric acid [m^2.kg.s^-2.mol^-1]
    aA_ref      = aA(muA_ref);                  // Reference activity of sulfuric acid
    aA_val      = aA(muA_val);                  // Activity of sulfuric acid
    i0_val      = i0(aA_ref, aA_val);           // Exchange current density [A.m^-2]
    Aa_val      = Aa(eps);                      // Specific interfacial area of lead [m^-1]
    sigmal_val  = sigmal(mA_val, eps);          // Effective ionic conductivity [s^3.A^2.m^-3.kg^-1]
    tH0_val     = tH0(mA_val);                  // Transference number of H+ with respect to the solvent velocity
    dphi0_val   = dphi0(muA_val, muA_ref);      // Interfacial voltage in equilibrium [m^2.kg.s^-3.A^-1]
    iBV_val     = iBV(i0_val, dphi0_val, dphi); // Butler-Volmer equation [A.m^-2]
    Dl_val      = Dl(eps, cA);                  // Effective diffusion coefficient of electrolyte [m^2.s^-1]
    V0_val      = V0(mA_val);                   // Partial molar volume of water [m^3.mol^-1]
    VA_val      = VA(mA_val);                   // Partial molar volume of sulfuric acid [m^3.mol^-1]



//-------------- Conservation of charge (interfacial voltage) -----------------
    
    // Flux1 = il
    Flux1       = zeros(N+1,1);
    Flux1(2:N)  = -0.5*(sigmal_val(1:N-1) + sigmal_val(2:N)).*(-(dphi(2:N) ..
                  - dphi(1:N-1)) + 0.5*(tH0_val(1:N-1) ..
                  + tH0_val(2:N) - 1).*(muA_val(2:N) - muA_val(1:N-1))/F)/dx;
    // Boundary condition (separator)
    Flux1(N+1)  = -Current(t)/S;
    // dF1 = grad(il)
    dF1         = (Flux1(2:N+1) - Flux1(1:N));
    
    // PDE 1
    r(1:N)      = -Aa_val*Cdl.*dphit + dF1/dx - Aa_val.*iBV_val;
    
    
//-------------- Conservation of volume (electrode porosity) ------------------
    
    // PDE 2
    r(N+1:N2)   = epst - 0.5/F*(MPb/rhoPb-MPbSO4/rhoPbSO4)*dF1/dx;
    
    
//------------- Conservation of species (acid concentration) ------------------
    
    // Flux2 = -Dl*cAx
    Flux2       = zeros(N+1,1);
    Flux2(2:N)  = -0.5*(Dl_val(1:N-1) + Dl_val(2:N)).*(cA(2:N) - cA(1:N-1))/dx;
    // Dirichlet right boundary (separator)
    cAr         = cA0 + Capacity(t)/(F*Vl);
    // Gradient of acid concentration [mol.m^-3]
    cAx         = zeros(N,1);
    cAx(2:N-1)  = 0.5*(cA(3:N) - cA(1:N-2))/dx;
    cANm12      = 0.5*(cA(N-1) + cA(N));    
    cAx(N)      = (cANm12 - 4*cA(N) + 3*cAr)/dx;
    // Boundary condition (separator)
    Flux2(N+1)  = -0.5*Dl_val(N)*cAx(N);
    // dF2 = grad(-Dl*cAx)
    dF2         = (Flux2(2:N+1) - Flux2(1:N));
    
    // Gradient of partial molar volume of sulfuric acid [m^3.mol^-1]
    VAx         = zeros(N,1);
    VAx(1)      = 0.5*(-VA_val(3) + 4*VA_val(2) - 3*VA_val(1))/dx;
    VAx(2:N-1)  = 0.5*(VA_val(3:N) - VA_val(1:N-2))/dx;
    VAx(N)      = 0.5*(VA_val(N-2) - 4*VA_val(N-1) + 3*VA_val(N))/dx;
    
    // Gradient of transference number of H+ with respect to the solvent velocity
    tH0x        = zeros(N,1);
    tH0x(1)     = 0.5*(-tH0_val(3) + 4*tH0_val(2) - 3*tH0_val(1))/dx;
    tH0x(2:N-1) = 0.5*(tH0_val(3:N) - tH0_val(1:N-2))/dx;
    tH0x(N)     = 0.5*(tH0_val(N-2) - 4*tH0_val(N-1) + 3*tH0_val(N))/dx;
    
    // PDE 3
    r(N2+1:N3)  = eps.*cAt + dF2/dx + (1 - cA.*VA_val)/F.*(Flux1(1:N).*tH0x) ..
                  + 0.5*(1 - cA.*VA_val).*(2*tH0_val - 1).*dF1/dx/F ..
                  - cA.*Dl_val./(1.0 - cA.*VA_val).*cAx.*VAx;
    
    // Without negligible terms
//    r(N2+1:N3)  = eps.*cAt + dF2/dx ..
//                  + 0.5*(1 - cA.*VA_val).*(2*tH0_val - 1).*dF1/dx/F;

    ires =  0;
end
