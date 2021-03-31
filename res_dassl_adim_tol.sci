function [r_, ires_] = res_dassl_adim_tol (t_, w_, wd_)
    
// --------------------------------------------------------- //
// RHS function                                              //
// Last modified: 02/04/2021                                 //
// --------------------------------------------------------- //

    global Call             // Declaration Call (for the calculation of work)
    Call = Call + 1;        // Count the number of function calls
    
    r_       = zeros(N3_,1);    // Initialization
    dphi_    = w_(1:N_);        // Dimensionless electrode voltage
    eps_     = w_(N_+1:N2_);    // Dimensionless porosity
    cA_      = w_(N2_+1:N3_);   // Dimensionless acid concentration
    dphit_   = wd_(1:N_);       // Dimensionless electrode voltage derived as a function of time
    epst_    = wd_(N_+1:N2_);   // Dimensionless porosity derived as a function of time
    cAt_     = wd_(N2_+1:N3_);  // Dimensionless acid concentration derived as a function of time


//------------------------ Physical relations ---------------------------------

    mA_val_      = mA_(cA_);                            // Dimensionless sulfuric acid molality
    muA_ref_     = muA_(mAref_);                        // Dimensionless reference electrochemical potential of sulfuric acid
    muA_val_     = muA_(mA_val_);                       // Dimensionless electrochemical potential of sulfuric acid
    aA_ref_      = aA_(muA_ref_);                       // Dimensionless reference activity of sulfuric acid
    aA_val_      = aA_(muA_val_);                       // Dimensionless activity of sulfuric acid
    i0_val_      = i0_(aA_ref_, aA_val_);               // Dimensionless exchange current density
    Aa_val_      = Aa_(eps_);                           // Dimensionless specific interfacial area of lead
    sigmal_val_  = sigmal_(mA_val_, eps_);              // Dimensionless effective ionic conductivity
    tH0_val_     = tH0_(mA_val_);                       // Dimensionless transference number of H+ with respect to the solvent velocity
    dphi0_val_   = dphi0_(muA_val_, muA_ref_);          // Dimensionless interfacial voltage in equilibrium
    iBV_val_     = iBV_(i0_val_, dphi0_val_, dphi_);    // Dimensionless Butler-Volmer equation
    Dl_val_      = Dl_(eps_, cA_);                      // Dimensionless effective diffusion coefficient of electrolyte
    V0_val_      = V0_(mA_val_);                        // Dimensionless partial molar volume of water
    VA_val_      = VA_(mA_val_);                        // Dimensionless partial molar volume of sulfuric acid



//-------------- Conservation of charge (interfacial voltage) -----------------
    
    // Flux1_ = il_
    Flux1_        = zeros(N_+1,1);
    Flux1_(2:N_)  = -0.5*(sigmal_val_(1:N_-1) + sigmal_val_(2:N_)).*(-(dphi_(2:N_) ..
                    - dphi_(1:N_-1)) + 0.5*(tH0_val_(1:N_-1) ..
                    + tH0_val_(2:N_) - 1).*(muA_val_(2:N_) - muA_val_(1:N_-1)))/dx_;
    // Boundary condition (separator)
    Flux1_(N_+1)  = -Current_(t_)/S_; // sep
    // dF1_ = grad_(il_)
    dF1_          = (Flux1_(2:N_+1) - Flux1_(1:N_)); // Flux1 = il
    
    // Dimensionless PDE 1
    r_(1:N_)      = -Aa_val_*Cdl_.*dphit_ + dF1_/dx_ - Aa_val_.*iBV_val_;
    
    
//-------------- Conservation of volume (electrode porosity) ------------------
    
    // Dimensionless PDE 2
    r_(N_+1:N2_)  = epst_ - 0.5*(MPb_/rhoPb_-MPbSO4_/rhoPbSO4_)*dF1_/dx_;
    
    
//------------- Conservation of species (acid concentration) ------------------
    
    // Flux2_ = -Dl_*cAx_
    Flux2_        = zeros(N_+1,1);
    Flux2_(2:N_)  = -0.5*(Dl_val_(1:N_-1) + Dl_val_(2:N_)).*(cA_(2:N_) - cA_(1:N_-1))/dx_;
    // Dirichlet right boundary (separator)
    cAr_          = cA0_ + Capacity_(t_)/Vl_;
    // Dimensionless gradient of acid concentration
    cAx_          = zeros(N_,1);
    cAx_(2:N_-1)  = 0.5*(cA_(3:N_) - cA_(1:N_-2))/dx_;
    cANm12_       = 0.5*(cA_(N_-1) + cA_(N_));
    cAx_(N_)      = (cANm12_ - 4*cA_(N_) + 3*cAr_)/dx_;
    // Boundary condition (separator)
    Flux2_(N_+1)  = -0.5*Dl_val_(N_)*cAx_(N_);
    // dF2_ = grad_(-Dl_*cAx_)
    dF2_          = (Flux2_(2:N_+1) - Flux2_(1:N_));
    
    // Dimensionless gradient of partial molar volume of sulfuric acid
    VAx_          = zeros(N_,1);
    VAx_(1)       = 0.5*(-VA_val_(3) + 4*VA_val_(2) - 3*VA_val_(1))/dx_;
    VAx_(2:N_-1)  = 0.5*(VA_val_(3:N_) - VA_val_(1:N_-2))/dx_;
    VAx_(N_)      = 0.5*(VA_val_(N_-2) - 4*VA_val_(N_-1) + 3*VA_val_(N_))/dx_;
    
    // Dimensionless gradient of transference number of H+ with respect to the solvent velocity
    tH0x_         = zeros(N_,1);
    tH0x_(1)      = 0.5*(-tH0_val_(3) + 4*tH0_val_(2) - 3*tH0_val_(1))/dx_;
    tH0x_(2:N_-1) = 0.5*(tH0_val_(3:N_) - tH0_val_(1:N_-2))/dx_;
    tH0x_(N_)     = 0.5*(tH0_val_(N_-2) - 4*tH0_val_(N_-1) + 3*tH0_val_(N_))/dx_;
    
    // Dimensionless PDE 3
    r_(N2_+1:N3_) = eps_.*cAt_ + dF2_/dx_ + (1 - cA_.*VA_val_).*(Flux1_(1:N_).*tH0x_) ..
                    + 0.5*(1 - cA_.*VA_val_).*(2*tH0_val_ - 1).*dF1_/dx_ ..
                    - cA_.*Dl_val_./(1.0 - cA_.*VA_val_).*cAx_.*VAx_;

    ires_ =  0;
end
