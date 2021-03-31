//< one line to give the library's name and a brief idea of what it does. >
//Copyright (C) < year > < name of author >
//This library is free software; 
//you can redistribute it and/or modify it under the terms of 
//the GNU Lesser General Public License as published by the Free Software Foundation; 
//either version 2.1 of the License, or (at your option) any later version.


// --------------------------------------------------------- //
// Definition of parameters and numerical solution in 1D     //
// battery model: one electrode configuration                //
// --------------------------------------------------------- //
// Last modified: 02/04/2021                                 //
// --------------------------------------------------------- //


clear;              // Kills variables
xdel(winsid());     // Clear all figures
clc;                // Clear Command Window
funcprot(0);        // No warning message when a function is redefined


// Load functions : Enter the path (no accent)
path = input("Enter the path (Ex: C:\Lead_Electrode) : ", "string");
// execution of the path
exec(path+'\res_dassl_adim.sci');

// Experimental data
load(path+"\data.sod");




//                          ////////////////////                               
//                          // Initialization //                               
//                          ////////////////////                               

//---------------------- Definition of parameters -----------------------------

alpha       = 0.5;          // Symmetry factor
epsc        = 0;            // Closed porosity
mu0Ag2SO4   = -618410.0;    // Gibbs free energy of formation of silver sulfate [m^2.kg.s^-2.mol^-1]
mu0H2SO4    = -690003.0;    // Gibbs free energy of formation of sulfuric acid [m^2.kg.s^-2.mol^-1]
mu0PbSO4    = -813140.0;    // Gibbs free energy of formation of lead-sulfate [m^2.kg.s^-2.mol^-1]
rhoAC       = 2380.0;       // Density of YP-50F activated carbon [kg.m^-3]
rhoDIA      = 2196.0;       // Density of diatomaceous earth (vitreous SiO2) [kg.m^-3]
rhoMCF      = 1800.0;       // Density of milled carbon fiber [kg.m^-3]
rhoPb       = 11341.0;      // Pycnometric density of pure lead [kg.m^-3]
rhoPbSO4    = 6291.0;       // Pycnometric density of lead-sulfate [kg.m^-3]
rhoVSA      = 500.0;        // Apparent density of sodium-lignosulfonate [kg.m^-3]
tau         = 1.03;         // Tortuosity factor
tau1        = 1.0 + tau;
Aa0         = 2.5e6;        // Specific interfacial area [m^-1]
cA0         = 4.12e3;       // Initial acid concentration [mol.m^-3]
Cdl         = 0.1;          // Double layer capacitance [s^4.A^2.m^-4.kg^-1]
F           = 96485.0;      // Faraday's constant [s.A.mol^-1]
k0          = 150;          // Exchange reaction rate [mol.s^-1.m^-2]
L           = 2.8e-4;       // Active material thickness [m]
mAref       = 6;            // Reference electrode sulfuric acid density of 1.26 kg/L
MPb         = 0.2072;       // Molecular weight of pure lead [kg.mol^-1]
MPbO        = 0.2232;       // Molecular weight of lead-monoxide [kg.mol^-1]
MPbSO4      = 0.3033;       // Molecular weight of lead-sulfate [kg.mol^-1]
R           = 8.3145;       // Universal gas constant [m^2.kg.s^-2.K^-1.mol^-1]
Rc          = 0;            // Connector resistance [m^2.kg.s^-3.A^-2]
S           = 2.75e-3;      // Electrode surface [m^2]
T0          = 295.0;        // Initial temperature [K]
T           = T0;           // Absolute temperature [K]
T2          = T*T;
Vl          = 1e-5;         // Volume of electrolyte [m^3]
Wam         = 3.051e-3;     // Mass of dried active material before formation [kg]
wAC         = 8.75e-3;      // Weight fraction of YP-50F activated carbon
wDIA        = 2.62e-2;      // Weight fraction of diatomaceous earth
wMCF        = 8.75e-3;      // Weight fraction of milled carbon fiber
wPbO        = 0.875;        // Weight fraction of lead-monoxide
wVSA        = 5.25e-3;      // Weight fraction of sodium-lignosulfonate



///////////// Definition of dimensionless parameters: ///////////////////////

// Reference parameters
//L          = 2.8e-4;
//T0         = 295.0;
//Wam        = 3.051e-3; 
//cA0        = 4.12e3;
//F          = 96485.0;
//R          = 8.3145;

// Dimensionless parameters
alpha_       = alpha;                               // Dimensionless symmetry factor
epsc_        = epsc;                                // Dimensionless closed porosity
mu0Ag2SO4_   = mu0Ag2SO4/(R*T0);                    // Dimensionless gibbs free energy of formation of silver sulfate
mu0H2SO4_    = mu0H2SO4/(R*T0);                     // Dimensionless gibbs free energy of formation of sulfuric acid
mu0PbSO4_    = mu0PbSO4/(R*T0);                     // Dimensionless gibbs free energy of formation of lead-sulfate
rhoAC_       = rhoAC*(L^3/Wam);                     // Dimensionless density of YP-50F activated carbon
rhoDIA_      = rhoDIA*(L^3/Wam);                    // Dimensionless density of diatomaceous earth (vitreous SiO2)
rhoMCF_      = rhoMCF*(L^3/Wam);                    // Dimensionless density of milled carbon fiber
rhoPb_       = rhoPb*(L^3/Wam);                     // Dimensionless pycnometric density of pure lead
rhoPbSO4_    = rhoPbSO4*(L^3/Wam);                  // Dimensionless pycnometric density of lead-sulfate
rhoVSA_      = rhoVSA*(L^3/Wam);                    // Dimensionless apparent density of sodium-lignosulfonate
tau_         = tau;                                 // Dimensionless tortuosity factor
tau1_        = 1.0 + tau_;
Aa0_         = Aa0*L;                               // Dimensionless specific interfacial area
cA0_         = cA0/cA0;                             // Dimensionless initial acid concentration
Cdl_         = Cdl*(R*T0/(cA0*L*F^2));              // Dimensionless double layer capacitance
k0_          = k0*sqrt(Wam/(R*T0*cA0^3*L^3));       // Dimensionless exchange reaction rate
mAref_       = mAref*(Wam/(cA0*L^3));               // Dimensionless reference electrode sulfuric acid density of 1.26 kg/L
MPb_         = MPb*(cA0*L^3/Wam);                   // Dimensionless molecular weight of pure lead
MPbO_        = MPbO*(cA0*L^3/Wam);                  // Dimensionless molecular weight of lead-monoxide
MPbSO4_      = MPbSO4*(cA0*L^3/Wam);                // Dimensionless molecular weight of lead-sulfate
Rc_          = Rc*(F^2*sqrt(cA0^3*L^7/(R*T0*Wam))); // Dimensionless connector resistance
S_           = S/L^2;                               // Dimensionless electrode surface
T_           = T/T0;                                // Dimensionless absolute temperature
T2_          = T_*T_;
Vl_          = Vl/L^3;                              // Dimensionless volume of electrolyte
wAC_         = wAC;                                 // Dimensionless weight fraction of YP-50F activated carbon
wDIA_        = wDIA;                                // Dimensionless weight fraction of diatomaceous earth
wMCF_        = wMCF;                                // Dimensionless weight fraction of milled carbon fiber
wPbO_        = wPbO;                                // Dimensionless weight fraction of lead-monoxide
wVSA_        = wVSA;                                // Dimensionless weight fraction of sodium-lignosulfonate


//-------------------- Dimensionless physical relations -----------------------

// Dimensionless exchange current density
function out = i0_(aAref_, aA_)
    out = 2.0*k0_*aAref_.^(1 - alpha_).*aA_.^alpha_;
endfunction

// Dimensionless sulfuric acid molality
function out = mA_(cA_)
    out = (Wam/L^3)*cA_.*(7.76e-4 + 4.77e-7*T0*T_ + cA0*cA_.*(4.48e-8 + 8.78e-11*T0*T_));
endfunction

// Dimensionless electrochemical potential of sulfuric acid
function out = muA_(mA_)
    out = (R*T0*mu0H2SO4_ + 12900.0*sqrt(cA0*L^3/Wam)*sqrt(mA_) - 79200.0)/(R*T0);
endfunction

// Dimensionless activity of sulfuric acid
function out = aA_(muA_)
    out = exp((muA_-mu0H2SO4_)/T_);
endfunction

// Dimensionless specific interfacial area of lead
function out = Aa_(eps_)

    out = Aa0_*((eps_-epsmin_)/(epsmax_-epsmin_)).^(2/3);
endfunction

// Dimensionless Butler-Volmer equation
function out = iBV_(i0_, dphi0_, dphi_)
    out = i0_.*(exp((1 - alpha_)*2.0/T_*(dphi_ - dphi0_)) - exp(-alpha_*2.0/T_*(dphi_ - dphi0_)));
endfunction

// Dimensionless interfacial voltage in equilibrium
function out = dphi0_(muA_, muAref_)
    out = 0.5*(mu0PbSO4_ - mu0Ag2SO4_ - muA_ + muAref_);
endfunction

// Dimensionless effective ionic conductivity
function out = sigmal_(mA_, eps_)
    out = sqrt(R*T0*Wam^3/(cA0^5*L^11*F^4))*(0.287*T0^2*T2_ - 119.0*T0*T_ + 1.25e4).*(eps_.^(tau1_)).*exp(-0.806*(log((cA0*L^3/Wam)*mA_) - 6.15e-3*T0*T_ - 1.2).^2)./mA_;
endfunction

// Dimensionless transference number of H+ with respect to the solvent velocity
function out = tH0_(mA_)
    out = (8.23e-6*T0*T_ - 2.71e-3)*(cA0^2*L^6/Wam^2).*mA_.^2 - (1.01e-4*T0*T_ - 1.61e-2)*(cA0*L^3/Wam).*mA_ - 9.88e-4*T0*T_ + 1.11;
endfunction

// Dimensionless effective diffusion coefficient of electrolyte
function out = Dl_(eps_, cA_)
    out = sqrt(Wam/(R*T0*cA0*L^5))*eps_.^(tau1_).*exp(-2174.0./(T0*T_)).*(2.65e-6 + 1.72e-2./(cA0*cA_).*exp(-0.587*(log(cA0*cA_) - 9.67).^2));
endfunction

// Dimensionless partial molar volume of water
function out = V0_(mA_)
    out = cA0*(-6.32e-8*(cA0*L^3/Wam)*mA_ + 1.81e-5);
endfunction

// Dimensionless partial molar volume of sulfuric acid
function out = VA_(mA_)
    out = cA0*(7.47e-7*(cA0*L^3/Wam)*mA_ + 3.92e-5);
endfunction


//----------------- Dimensionless initial interfacial voltage -----------------

// Dimensionless initial Sulfuric acid molality
mA0_     = mA_(cA0_);
// Dimensionless refrence electrochemical potential of sulfuric acid
muA_ref_ = muA_(mAref_);
// Dimensionless initial electrochemical potential of sulfuric acid
muA0_    = muA_(mA0_);
// Dimensionless initial interfacial voltage in equilibrium
dphi00_  = dphi0_(muA0_, muA_ref_);


//------------------------------ Porosity -------------------------------------

// Dimensionless volume fraction of inert materials
epsi_      = (wAC_/rhoAC_+wDIA_/rhoDIA_+wMCF_/rhoMCF_+wVSA_/rhoVSA_)/S_ + epsc_;
// Dimensionless initial volume fraction of lead
epsPb0_    = MPb_*wPbO_/(MPbO_*rhoPb_)/S_;
// Dimensionless initial volume fraction of lead sulfate
epsPbSO40_ = 0;
// Dimensionless initial electrode porosity (filled with electrolyte)
eps0_      = 1-epsPb0_-epsPbSO40_-epsi_;
// Dimensionless domain of definition of electrode porosity
epsmax_    = eps0_;
epsmin_    = 0.327; 


//------------------Dimensionless numerical parameters ------------------------

N_  = 10;                        // Number of discretization points
N2_ = 2*N_; 
N3_ = 3*N_;                     // Total of 3*N discretization points (3 equations to solve)
x_  = linspace(0, L, N_)'/L;    // Dimensionless space discretization (electrode thickness)
dx_ = x_(2) - x_(1);            // Dimensionless space variation (constant step)




//                         /////////////////                                    
//                         // SIMULATIONS //                                    
//                         /////////////////                                   

// Capacity regimes studied (20h discharge - 5min discharge)
rate     = ["0.05C","0.2C","0.5C","1C","2C","4C","8C","12C"];
rate_EXP = ["C20","C5","C2","1C","2C","4C","8C","12C"];
// Color map for display
cmap = flipdim([1,0,0;1,0.5,0;0.9,0.9,0;0,1,0;0,1/2,0;0,1,1;0,0,1;0.7,0,0.7],1);
rtol = 1e-6         // Relative tolerance
atol = 1e-7         // Absolute tolerance
ti   = getdate();   // Computation time initialization
mprintf("\nRate\tMode\tComput. time\tElapsed time\tRemaining time\n");
EXP  = list();      // Experimental data
SIM_ = list();      // Simulated data
figure('Figure_name','Uneg_wrt_Q','BackgroundColor',[1 1 1]);

// Simulation for each capacity regime
for ii = 1:size(rate,2)
//                       ///////////////////////                               
//                       // Experimental data //                               
//                       ///////////////////////                               

    // Recovery of experimental data
    execstr("EXP($+1) = DCH.i"+rate_EXP(ii)+";"); 
    // Dimensionless time and capacity (experimental data)
    tEXP_ = sqrt(R*T0*cA0*L/Wam)*EXP(ii).t;
    QEXP_ = EXP(ii).Q/(F*cA0*L^3);
    // Dimensionless polynomial interpolation of capacity
    function out = Capacity_(t_)
        out = interp(t_,tEXP_,QEXP_,splin(tEXP_,QEXP_),"natural");
    endfunction


//                      //////////////////////////                             
//                      // Discharge simulation //                             
//                      //////////////////////////                             


    // Discharge data
    EOD(ii) = find(EXP(ii).i == 0,2)($)-1;          // End of discharge
    t_DCH   = EXP(ii).t(1:EOD(ii));                 // DCH time (s)
    i_DCH   = EXP(ii).i(1:EOD(ii));                 // DCH current (A)
    t_DCH_  = t_DCH*sqrt(R*T0*cA0*L/Wam);           // Dimensionless DCH time
    i_DCH_  = i_DCH*sqrt(Wam/(R*T0*cA0^3*L^7*F^2)); // Dimensionless DCH current
    // Dimensionless polynomial interpolation of discharge current
    function out = Current_(t_)
        out = interp(t_,t_DCH_,i_DCH_,splin(t_DCH_,i_DCH_),"natural");
    endfunction
    // Computation
    w0_            = zeros(N3_,1);                  // Variable initial conditions
    w0_(1:N_)      = dphi00_;                       // Dimensionless voltage initial condition
    w0_(N_+1:N2_)  = eps0_;                         // Dimensionless porosity initial condition
    w0_(N2_+1:N3_) = cA0_;                          // Dimensionless concentration initial condition
    wd0_           = zeros(N3_,1);                  // Time derivative of var. initial conditions
    t0_            = t_DCH_(1);                     // Dimensionless initial time for the simulation
    t1             = getdate();                     // Start execution time calculation
    W_DCH_         = dassl([w0_, wd0_], t0_, t_DCH_,...
                           rtol, atol, res_dassl_adim);
    W_DCH_         = W_DCH_(2:$,1:$);               // 1st line: time
    t2             = getdate();                     // Finish execution time calculation

    // Display of execution time (discharge)
    mprintf(rate(ii)+"\tDCH\t%i min %.3f s\t%i min %.3f s\t%i min %i s\n",..
    int8(etime(t2,t1)/60),modulo(etime(t2,t1),60),..
    int8(etime(t2,ti)/60),modulo(etime(t2,ti),60),..
    int8((etime(t2,ti)/((ii-0.5)/size(rate,2))-etime(t2,ti))/60),..
    modulo(etime(t2,ti)/((ii-0.5)/size(rate,2))-etime(t2,ti),60));


//                        //////////////////////                               
//                        // Pause simulation //                               
//                        //////////////////////                               
    // Pause data
    t_PAU  = EXP(ii).t(EOD(ii):$);                  // PAU time (s)
    i_PAU  = EXP(ii).i(EOD(ii):$);                  // PAU current (A)
    // Dimensionless time and current
    t_PAU_ = t_PAU*sqrt(R*T0*cA0*L/Wam);            // Dimensionless PAU time
    i_PAU_ = i_PAU*sqrt(Wam/(R*T0*cA0^3*L^7*F^2));  // Dimensionless PAU current
    // Dimensionless polynomial interpolation of pause current
    function out = Current_(t_)
        out = interp(t_,t_PAU_,i_PAU_,splin(t_PAU_,i_PAU_),"natural");
    endfunction

    // Computation
    w0_    = W_DCH_(1:N3_,$);                       // Dimensionless solution at the last discharge time step
    wd0_   = W_DCH_(N3_+1:$,$);                     // Dimensionless derivative solution at the last discharge time step
    t0_    = t_PAU_(1);                             // Dimensionless initial time for the simulation
    t1     = getdate();                             // Start execution time calculation
    // Problem solving by DASSL solver (pause)
    W_PAU_ = dassl([w0_, wd0_], t0_, t_PAU_, ...
                   rtol, atol, res_dassl_adim);
    W_PAU_ = W_PAU_(2:$,1:$);                       // 1st line: time
    t2     = getdate();                             // Finish execution time calculation

    // Display of execution time (pause)
    mprintf(rate(ii)+"\tPAU\t%i min %.3f s\t%i min %.3f s\t%i min %i s\n",..
    int8(etime(t2,t1)/60),modulo(etime(t2,t1),60),..
    int8(etime(t2,ti)/60),modulo(etime(t2,ti),60),..
    int8((etime(t2,ti)/(ii/size(rate,2))-etime(t2,ti))/60),..
    modulo(etime(t2,ti)/(ii/size(rate,2))-etime(t2,ti),60));
    
    
    
//                   /////////////////////////////                            
//                   // Save simulation results //                            
//                   /////////////////////////////                            
    SIM_($+1).t    = EXP(ii).t;                                             // Dimensionless time
    SIM_(ii).Q     = EXP(ii).Q;                                             // Dimensionless capacity
    SIM_(ii).dphi  = [W_DCH_(1:N_,:)*(R*T0/F),W_PAU_(1:N_,2:$)*(R*T0/F)];   // Dimensionless simulated electrode voltage (DCH and PAU)
    SIM_(ii).eps   = [W_DCH_(N_+1:N2_,:),W_PAU_(N_+1:N2_,2:$)];             // Dimensionless simulated porosity (DCH and PAU)
    SIM_(ii).cA    = [W_DCH_(N2_+1:N3_,:)*cA0,W_PAU_(N2_+1:N3_,2:$)*cA0];   // Dimensionless simulated acid concentration (DCH and PAU)
    SIM_(ii).u_neg = SIM_(ii).dphi(1,:)' + Rc*EXP(ii).i;                    // Dimensionless electrode voltage



//                   /////////////////////////////                            
//                   // Save simulation results //                            
//                   /////////////////////////////                            
    // Electrode voltage as a function of capacity (plot for each regime)
    // Experimental data
    plot(-EXP(ii).Q(1:EOD(ii))/3600,EXP(ii).u_neg(1:EOD(ii)),"LineS","--",..
    "color",cmap(ii,:));
    // Simulated data
    plot(-SIM_(ii).Q(1:EOD(ii))/3600,SIM_(ii).u_neg(1:EOD(ii)),"color",cmap(ii,:),..
    "thick",2);
end

// Figure parameters
a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on"; a.font_size = 3.3;
// Axis
a.x_label.text = '$\text{Capacity, }Q\text{ (Ah)}$';
a.x_label.font_size = 3.3; a.x_ticks.labels = '$'+a.x_ticks.labels+'$';
a.y_label.text = '$\text{Electrode voltage, }u\text{ (V)}$';
a.y_label.font_size = 3.3; a.y_ticks.labels = '$'+a.y_ticks.labels+'$';
// Legend
hl = legend('$\text{0.05C exp}$','$\text{0.05C sim}$','$\text{0.2C exp}$',..
'$\text{0.2C sim}$','$\text{0.5C exp}$','$\text{0.5C sim}$','$\text{1C exp}$',..
'$\text{1C sim}$','$\text{2C exp}$','$\text{2C sim}$','$\text{4C exp}$',..
'$\text{4C sim}$','$\text{8C exp}$','$\text{8C sim}$','$\text{12C exp}$',..
'$\text{12C sim}$');
hl.thickness = 0; hl.background = 8; hl.line_mode = 'off'; hl.font_color = -1;
hl.legend_location = "out_upper_right"; hl.font_size = 3;


//                   /////////////////////////////                            
//                   // Plot simulation results //                            
//                   /////////////////////////////                            

// Display electrode voltage pause as a function of time
figure('Figure_name','Uneg_wrt_t_PAU','BackgroundColor',[1 1 1]);
f = gcf(); f.figure_size = [626*1.5,587*2];
// plot for each regime
for ii = 1:size(rate,2)
    EOD(ii) = find(EXP(ii).i == 0,2)($)-1;  // End of discharge
    t_PAU   = EXP(ii).t(EOD(ii):$);         // PAU time (s)
    // Experimental data
    plot((t_PAU-t_PAU(1))/60,EXP(ii).u_neg(EOD(ii):$),"LineS","--","color",cmap(ii,:));
    // Simulated data
    plot((t_PAU-t_PAU(1))/60,SIM_(ii).u_neg(EOD(ii):$),"color",cmap(ii,:),"thick",2);
end
// Figure parameters
a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on";
a.font_size = 3.3; a.data_bounds(:,1) = [0;30]; a.data_bounds(:,2) = [-1;-0.8];
// Axis
a.x_label.text = '$\text{Time, }t\text{ (min)}$';
a.x_ticks.labels = '$'+a.x_ticks.labels+'$'; a.x_label.font_size = 3.3;
a.y_label.text = '$\text{Voltage, }u\text{ (V)}$';
a.y_ticks.labels = '$'+a.y_ticks.labels+'$'; a.y_label.font_size = 3.3;
// Legend
hl = legend('$\text{0.05C exp}$','$\text{0.05C sim}$','$\text{0.2C exp}$',..
'$\text{0.2C sim}$','$\text{0.5C exp}$','$\text{0.5C sim}$','$\text{1C exp}$',..
'$\text{1C sim}$','$\text{2C exp}$','$\text{2C sim}$','$\text{4C exp}$',..
'$\text{4C sim}$','$\text{8C exp}$','$\text{8C sim}$','$\text{12C exp}$',..
'$\text{12C sim}$');
hl.thickness = 0; hl.background = 8; hl.line_mode = 'off'; hl.font_color = -1;
hl.legend_location = "in_upper_right"; hl.font_size = 3;


// Display electrode voltage as a function of time
figure('Figure_name','Uneg_wrt_t','BackgroundColor',[1 1 1]);
f = gcf(); f.figure_size = [626*1.5,587*2];
for ii = 1:size(rate,2)
    subplot(2,4,ii);
    if ii < 5 then
        plot(EXP(ii).t/3600,EXP(ii).u_neg,"LineS","--","color",cmap(ii,:));
        plot(SIM_(ii).t/3600,SIM_(ii).u_neg,"color",cmap(ii,:),"thick",2);
        a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on";
        a.font_size = 3.3; a.data_bounds(:,2) = [-1;-0.2];
        a.x_label.text = '$\text{Time, }t\text{ (h)}$';
    else
        plot(EXP(ii).t/60,EXP(ii).u_neg,"LineS","--","color",cmap(ii,:));
        plot(SIM_(ii).t/60,SIM_(ii).u_neg,"color",cmap(ii,:),"thick",2);
        a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on";
        a.font_size = 3.3; a.data_bounds(:,2) = [-1;-0.2];
        a.x_label.text = '$\text{Time, }t\text{ (min)}$';
    end
    a.x_ticks.labels = '$'+a.x_ticks.labels+'$'; a.x_label.font_size = 3.3;
    a.y_label.text = '$\text{Voltage, }u\text{ (V)}$';
    a.y_ticks.labels = '$'+a.y_ticks.labels+'$'; a.y_label.font_size = 3.3;
    hl = legend('$\text{'+rate(ii)+' exp}$','$\text{'+rate(ii)+' sim}$');
    hl.thickness = 0; hl.background = 8; hl.line_mode = 'off'; hl.font_color = -1;
    hl.legend_location = "in_upper_right"; hl.font_size = 3.3;
end


// Display electrode porosity according to its thickness
figure('Figure_name','eps_DCH','BackgroundColor',[1 1 1]);
f = gcf(); f.figure_size = [626*2,587];
for ii = 1:size(rate,2)
    subplot(2,4,ii);
    // colormap with a color for each time step
    cmap = jetcolormap(EOD(ii)); 
    for jj = 1:ceil(EOD(ii)/100):EOD(ii)
        plot(1e3*x_,SIM_(ii).eps(:,jj),"color",cmap(jj,:));
    end
    a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on";
    a.font_size = 3.3; a.data_bounds(:,2) = [0.3;0.65];
    a.x_label.text = '$\text{Electrode thickness, }x\text{ (mm)}$';
    a.x_ticks.labels = '$'+a.x_ticks.labels+'$'; a.x_label.font_size = 3.3;
    a.y_label.text = '$\text{Porosity, }\epsilon$';
    a.y_ticks.labels = '$'+a.y_ticks.labels+'$'; a.y_label.font_size = 3.3;
end


// Display acid concentration according to electrode thickness
figure('Figure_name','cA_DCH','BackgroundColor',[1 1 1]);
f = gcf(); f.figure_size = [626*2,587];
for ii = 1:size(rate,2)
    subplot(2,4,ii);
    // colormap with a color for each time step
    cmap = jetcolormap(EOD(ii)); 
    for jj = 1:ceil(EOD(ii)/100):EOD(ii)
        plot(1e3*x_,SIM_(ii).cA(:,jj),"color",cmap(jj,:));
    end
    a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on";
    a.font_size = 3.3; a.data_bounds(:,2) = [2500;4300];
    a.x_label.text = '$\text{Electrode thickness, }x\text{ (mm)}$';
    a.x_ticks.labels = '$'+a.x_ticks.labels+'$'; a.x_label.font_size = 3.3;
    a.y_label.text = '$\text{Concentration, }c_\mathrm{A}\text{ (mol.m^{-3})}$';
    a.y_ticks.labels = '$'+a.y_ticks.labels+'$'; a.y_label.font_size = 3.3;
end
