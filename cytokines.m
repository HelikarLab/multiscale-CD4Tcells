function [FinalC1, FinalC2, FinalC3]=cytokines(Ctref, Clnref, Ccref, Ct0, Cln0, Cc0, A1, A2, A3, Nt, Nln, Nc, tau, Cs, kdeg, Qa, Qb, Vt, Vln, Vc, Pt, Pt_in, Pt_ex, Pln, Pln_in, Pln_ex, Pc, Pc_in, Pc_ex)
    %This function calculates the concentrations of a cytokine at the end of
    %a time step. The inputs are all from the beginning of the time step.

    %Define the dimensionless parameters in the cytokine ODEs.
    L1=Pt*tau/Cs;
    L2=Pt_in*tau/Cs;
    L3=Pt_ex*tau/Cs;
    L4=kdeg*tau;
    L5=Qa*tau/Vt;
    L6=Pln*tau/Cs;
    L7=Pln_in*tau/Cs;
    L8=Pln_ex*tau/Cs;
    L9=Qb*tau/Vln;
    L10=Pc*tau/Cs;
    L11=Pc_in*tau/Cs;
    L12=Pc_ex*tau/Vc;
    L13=(Qa-Qb)*tau/Vc;
    L14=Qb*tau/Vc;
    L15=Qa*tau/Vc;
    
    %Define the cytokine ODEs, set the number of time discretisation steps,
    %and solve for the concentrations after one time step.
    RHS=@(t, C) [L1+L2*A1+L3*Nt-L4*C(1)+L5*(C(3)-C(1)); L6+L7*A2+L8*Nln-L4*C(2)+L9*(C(1)-C(2)); L10+L11*A3+L12*Nc-L4*C(3)+L13*C(1)+L14*C(2)-L15*C(3)];
    dis_steps=10;
    [t, C]=ode45(RHS, [0:1/dis_steps:1], [Ct0, Cln0, Cc0]);
    
    %Return the concentrations at the end of the time step if they are
    %above the baselines. Otherwise, return the baseline concentrations.
    if C(dis_steps+1, 1)>Ctref
        FinalC1=C(dis_steps+1, 1);
    else
        FinalC1=Ctref;
    end
    
    if C(dis_steps+1, 2)>Clnref
        FinalC2=C(dis_steps+1, 2);
    else
        FinalC2=Clnref;
    end
    
    if C(dis_steps+1, 3)>Ccref
        FinalC3=C(dis_steps+1, 3);
    else
        FinalC3=Ccref;
    end
end