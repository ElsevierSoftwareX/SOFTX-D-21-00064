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
global Call;        // Declaration Call (for the calculation of work)


// Load functions : Enter the path (no accent)
path = input("Enter the path (Ex: C:\Lead_Electrode) : ", "string");
// execution of the path
exec(path+'\res_dassl_dim_tol.sci');

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
mAref       = 6;            // Reference electrode sulfuric acid density of 1.26 kg/L [mol.kg^-1]
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


//------------------------- Physical relations --------------------------------

// Exchange current density [A.m^-2]
function out = i0(aAref, aA)
    out = 2.0*F*k0*aAref.^(1 - alpha).*aA.^alpha;
endfunction

// Sulfuric acid molality [mol.kg^-1]
function out = mA(cA)
    out = cA.*(7.76e-4 + 4.77e-7*T + cA.*(4.48e-8 + 8.78e-11*T));
endfunction

// Electrochemical potential of sulfuric acid [m^2.kg.s^-2.mol^-1]
function out = muA(mA)
    out = mu0H2SO4 + 12900.0*sqrt(mA) - 79200.0;
endfunction

// Activity of sulfuric acid
function out = aA(muA)
    out = exp((muA-mu0H2SO4)./(R*T));
endfunction

// Specific interfacial area of lead [m^-1]
function out = Aa(eps)
    out = Aa0*((eps-epsmin)/(epsmax-epsmin)).^(2/3);
endfunction

// Butler-Volmer equation [A.m^-2]
function out = iBV(i0, dphi0, dphi)
    out = i0.*(exp((1 - alpha)*2.0*F/(R*T)*(dphi - dphi0)) ...
               - exp(-alpha*2.0*F/(R*T)*(dphi - dphi0)));
endfunction

// Interfacial voltage in equilibrium [m^2.kg.s^-3.A^-1]
function out = dphi0(muA, muAref)
    out = 0.5*(mu0PbSO4 - mu0Ag2SO4 - muA + muAref)/F;
endfunction

// Effective ionic conductivity [s^3.A^2.m^-3.kg^-1]
function out = sigmal(mA, eps)
    out = (0.287*T2 - 119.0*T + 1.25e4).*(eps.^(tau1)) ...
           .*exp(-0.806*(log(mA) - 6.15e-3*T - 1.2).^2)./mA;
endfunction

// Transference number of H+ with respect to the solvent velocity
function out = tH0(mA)
    out = (8.23e-6*T - 2.71e-3).*mA.^2 - (1.01e-4*T - 1.61e-2).*mA ...
           - 9.88e-4*T + 1.11;
endfunction

// Effective diffusion coefficient of electrolyte [m^2.s^-1]
function out = Dl(eps, cA)
    out = eps.^(tau1).*exp(-2174.0./T).*(2.65e-6 + (1.72e-2./cA) ...
          .*exp(-0.587*(log(cA) - 9.67).^2));
endfunction

// Partial molar volume of water [m^3.mol^-1]
function out = V0(mA)
    out = -6.32e-8*mA + 1.81e-5;
endfunction

// Partial molar volume of sulfuric acid [m^3.mol^-1]
function out = VA(mA)
    out = 7.47e-7*mA + 3.92e-5;
endfunction


//----------------------- Initial interfacial voltage -------------------------

// Initial Sulfuric acid molality [mol.kg^-1]
mA0     = mA(cA0);
// Refrence electrochemical potential of sulfuric acid [m^2.kg.s^-2.mol^-1]
muA_ref = muA(mAref);
// Initial electrochemical potential of sulfuric acid [m^2.kg.s^-2.mol^-1]
muA0    = muA(mA0);
// Initial interfacial voltage in equilibrium [m^2.kg.s^-3.A^-1]
dphi00  = dphi0(muA0, muA_ref);


//------------------------------ Porosity -------------------------------------

// Volume fraction of inert materials
epsi      = Wam/(L*S)*(wAC/rhoAC+wDIA/rhoDIA+wMCF/rhoMCF+wVSA/rhoVSA) + epsc;
// Initial volume fraction of lead
epsPb0    = Wam*MPb*wPbO/(L*S*MPbO*rhoPb);
// Initial volume fraction of lead sulfate
epsPbSO40 = 0;
// Initial electrode porosity (filled with electrolyte)
eps0      = 1-epsPb0-epsPbSO40-epsi;
// Domain of definition of electrode porosity
epsmax    = eps0;
epsmin    = 0.327;


//------------------------ Numerical parameters -------------------------------

N  = 10;                    // Number of discretization points
N2 = 2*N;
N3 = 3*N;                   // Total of 3*N discretization points (3 equations to solve)
x  = linspace(0, L, N)';    // Space discretization (electrode thickness)
dx = x(2) - x(1);           // Space variation (constant step)




//                         /////////////////                                    
//                         // SIMULATIONS //                                    
//                         /////////////////                                   

// Capacity regimes studied (20h discharge - 5min discharge)
rate     = ["0.05C","0.2C","0.5C","1C","2C","4C","8C","12C"];
rate_EXP = ["C20","C5","C2","1C","2C","4C","8C","12C"];
// Color map for display
cmap   = flipdim([1,0,0;1,0.5,0;0.9,0.9,0;0,1,0;0,1/2,0;0,1,1;0,0,1;0.7,0,0.7],1);
v_rtol = [1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11];    // Vector of relative tolerance
v_atol = [1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12];   // Vector of absolute tolerance
EXP    = list();                                    // Experimental data
SIM    = list();                                    // Simulated data
figure('Figure_name','Uneg_wrt_Q','BackgroundColor',[1 1 1]);

// Simulation for each tolerance
for idx_tol = 1:length(v_rtol)
    rtol = v_rtol(idx_tol);             // Relative tolerance
    atol = v_atol(idx_tol);             // Absolute tolerance
    ti     = getdate();                 // Computation time initialization
    mprintf("\nRelative tolerance: %.1e",rtol);
    mprintf("\nRate\tMode\tComput. time\tElapsed time\tRemaining time\n");
    
    // Simulation for each capacity regime
    for ii = 1:size(rate,2)
//                       ///////////////////////                               
//                       // Experimental data //                               
//                       ///////////////////////                               

        // Recovery of experimental data
        execstr("EXP($+1) = DCH.i"+rate_EXP(ii)+";");
        // Polynomial interpolation of capacity
        function out = Capacity(t)
            out = interp(t,EXP(ii).t,EXP(ii).Q,splin(EXP(ii).t,EXP(ii).Q),"natural");
        endfunction
        
//                      //////////////////////////                             
//                      // Discharge simulation //                             
//                      //////////////////////////                             

        Call = 0;                               // Resetting the work variable (DCH)
        // Discharge data
        EOD(ii) = find(EXP(ii).i == 0,2)($)-1;  // End of discharge
        t_DCH   = EXP(ii).t(1:EOD(ii));         // DCH time (s)
        i_DCH   = EXP(ii).i(1:EOD(ii));         // DCH current (A)
        // Polynomial interpolation of discharge current
        function out = Current(t)
            out = interp(t,t_DCH,i_DCH,splin(t_DCH,i_DCH),"natural");
        endfunction
        // Computation
        w0          = zeros(N3,1);              // Variable initial conditions
        w0(1:N)     = dphi00;                   // Voltage initial condition
        w0(N+1:N2)  = eps0;                     // Porosity initial condition
        w0(N2+1:N3) = cA0;                      // Concentration initial condition
        wd0         = zeros(N3,1);              // Time derivative of var. initial conditions
        t0          = t_DCH(1);                 // Initial time for the simulation
        t1          = getdate();                // Start execution time calculation
        // Problem solving by DASSL solver (discharge)
        W_DCH       = dassl([w0, wd0], t0, t_DCH, ...
                            rtol, atol, res_dassl_dim_tol);
        W_DCH       = W_DCH(2:$,1:$);           // 1st line: time
        t2          = getdate();                // Finish execution time calculation
        Call_DCH(ii,idx_tol) = Call;            // Total number of RHS function calls for each tolerance (DCH)
        
        // Display of execution time (discharge)
        mprintf(rate(ii)+"\tDCH\t%i min %.2f s\t%i min %.2f s\t%i min %i s\n",..
        int8(etime(t2,t1)/60),modulo(etime(t2,t1),60),..
        int8(etime(t2,ti)/60),modulo(etime(t2,ti),60),..
        int8((etime(t2,ti)/((ii-0.5)/size(rate,2))-etime(t2,ti))/60),..
        modulo(etime(t2,ti)/((ii-0.5)/size(rate,2))-etime(t2,ti),60));



//                        //////////////////////                               
//                        // Pause simulation //                               
//                        //////////////////////                               

        Call  = 0;                          // Resetting the work variable (PAU)
        // Pause data
        t_PAU = EXP(ii).t(EOD(ii):$);       // PAU time (s)
        i_PAU = EXP(ii).i(EOD(ii):$);       // PAU current (A)
        // Polynomial interpolation of pause current
        function out = Current(t)
            out = interp(t,t_PAU,i_PAU,splin(t_PAU,i_PAU),"natural");
        endfunction
        // Computation
        w0    = W_DCH(1:N3,$);              // Solution at the last discharge time step
        wd0   = W_DCH(N3+1:$,$);            // Derivative solution at the last discharge time step
        t0    = t_PAU(1);                   // Initial time for the simulation
        t1    = getdate();                  // Start execution time calculation
        // Problem solving by DASSL solver (pause)
        W_PAU = dassl([w0, wd0], t0, t_PAU, rtol, atol, res_dassl_dim_tol);
        W_PAU = W_PAU(2:$,1:$);             // 1st line: time
        t2    = getdate();                  // Finish execution time calculation
        Call_PAU(ii,idx_tol) = Call;        // Total number of RHS function calls for each tolerance (PAU)

        // Display of execution time (pause)
        mprintf(rate(ii)+"\tPAU\t%i min %.2f s\t%i min %.2f s\t%i min %i s\n",..
        int8(etime(t2,t1)/60),modulo(etime(t2,t1),60),..
        int8(etime(t2,ti)/60),modulo(etime(t2,ti),60),..
        int8((etime(t2,ti)/(ii/size(rate,2))-etime(t2,ti))/60),..
        modulo(etime(t2,ti)/(ii/size(rate,2))-etime(t2,ti),60));
    
    
    
//                   /////////////////////////////                            
//                   // Save simulation results //                            
//                   /////////////////////////////                            
        SIM($+1).t    = EXP(ii).t;                              // Time
        SIM(ii).Q     = EXP(ii).Q;                              // Capacity
        SIM(ii).dphi  = [W_DCH(1:N,:),W_PAU(1:N,2:$)];          // Simulated electrode voltage (DCH and PAU)
        SIM(ii).eps   = [W_DCH(N+1:N2,:),W_PAU(N+1:N2,2:$)];    // Simulated porosity (DCH and PAU)
        SIM(ii).cA    = [W_DCH(N2+1:N3,:),W_PAU(N2+1:N3,2:$)];  // Simulated acid concentration (DCH and PAU)
        SIM(ii).u_neg = SIM(ii).dphi(1,:)' + Rc*EXP(ii).i;      // Electrode voltage 



//                   /////////////////////////////                            
//                   // Plot simulation results //                            
//                   /////////////////////////////                            

        // Electrode voltage as a function of capacity (plot for each regime)
        // plot only for the first tolerance
        if idx_tol == 1 then
            // Experimental data
            plot(-EXP(ii).Q(1:EOD(ii))/3600,EXP(ii).u_neg(1:EOD(ii)),"LineS","--",..
            "color",cmap(ii,:));
            // Simulated data
            plot(-SIM(ii).Q(1:EOD(ii))/3600,SIM(ii).u_neg(1:EOD(ii)),"color",cmap(ii,:),..
            "thick",2);
        end
    end
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
    plot((t_PAU-t_PAU(1))/60,SIM(ii).u_neg(EOD(ii):$),"color",cmap(ii,:),"thick",2);
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
        plot(SIM(ii).t/3600,SIM(ii).u_neg,"color",cmap(ii,:),"thick",2);
        a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on";
        a.font_size = 3.3; a.data_bounds(:,2) = [-1;-0.2];
        a.x_label.text = '$\text{Time, }t\text{ (h)}$';
    else
        plot(EXP(ii).t/60,EXP(ii).u_neg,"LineS","--","color",cmap(ii,:));
        plot(SIM(ii).t/60,SIM(ii).u_neg,"color",cmap(ii,:),"thick",2);
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
        plot(1e3*x,SIM(ii).eps(:,jj),"color",cmap(jj,:));
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
        plot(1e3*x,SIM(ii).cA(:,jj),"color",cmap(jj,:));
    end
    a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on";
    a.font_size = 3.3; a.data_bounds(:,2) = [2500;4300];
    a.x_label.text = '$\text{Electrode thickness, }x\text{ (mm)}$';
    a.x_ticks.labels = '$'+a.x_ticks.labels+'$'; a.x_label.font_size = 3.3;
    a.y_label.text = '$\text{Concentration, }c_\mathrm{A}\text{ (mol.m^{-3})}$';
    a.y_ticks.labels = '$'+a.y_ticks.labels+'$'; a.y_label.font_size = 3.3;
end


// Color map for display
cmap = flipdim([1,0,0;1,0.5,0;0.9,0.9,0;0,1,0;0,1/2,0;0,1,1;0,0,1;0.7,0,0.7],1);
// Display Number of iterations in terms of tolerance
Marker = ["s","d","+","o","*","x","p","v"];     // varies the marker
figure('Figure_name','Work - Tolerance','BackgroundColor',[1 1 1]);
// Plot for each regime
for ii = 1:size(rate,2)
    // Discharge state
    plot("ll",v_rtol,Call_DCH(ii,:),'marker',Marker(ii),"LineS","-","color",cmap(ii,:));
end
for ii = 1:size(rate,2)
    // Pause state
    plot("ll",v_rtol,Call_PAU(ii,:),'marker',Marker(ii),"LineS","--","color",cmap(ii,:));
end
// Figure parameters
a = gca(); col=color("grey"); a.grid = [col,col]; a.box = "on"; a.font_size = 3.3;
a.data_bounds(:,2) = [1e2;1e5]; a.sub_ticks = [0,0];
// Axis
a.x_label.text = '$\text{Relative Tolerance}$';
a.x_label.font_size = 3.3; a.x_ticks.labels = '$10^{'+string(evstr(part(a.x_ticks.labels,3:5)))+'}$';
a.y_label.text = '$\text{Number of function calls}$';
a.y_label.font_size = 3.3; a.y_ticks.labels = '$10^{'+string(log10(evstr(strsubst(a.y_ticks.labels," ",""))))+'}$';
// Legend
hl = legend('$\text{C/20}$','$\text{C/5}$','$\text{C/2}$','$\text{1C}$',..
            '$\text{2C}$','$\text{4C}$','$\text{8C}$','$\text{12C}$',...
            '$\text{C/20}$','$\text{C/5}$','$\text{C/2}$','$\text{1C}$',..
            '$\text{2C}$','$\text{4C}$','$\text{8C}$','$\text{12C}$');
hl.thickness = 0; hl.background = 8; hl.line_mode = 'off'; hl.font_color = -1;
hl.legend_location = "out_upper_right"; hl.font_size = 3;
