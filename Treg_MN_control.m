function output = Treg_MN_control(Treg_model, Phase, solution, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, BN_glycolysis, BN_glu_uptake, BN_aa_tran, BN_mit_ox, BN_lip_eff, BN_glutaminolysis, BN_lip_syn)
    %This function uses the outputs of the BN model to parametrise the MN
    %model for phenotype 15 (Treg) and then calculates various production
    %rates.
    
    %Make a copy of the metabolic model for Treg.
    model1=Treg_model;
    
    %Use the converted outputs from the BN model to control the upper and lower
    %bounds of different classes of metabolic fluxes by knocking down the
    %optimised fluxes. For example, if the converted BN node for glycolysis (lactase)
    %is 0.2 at a time step, an optimised flux relevant to the process
    %is 100 and the bounds are -1000 and 1000, set its upper bound to 20
    %and the lower bound to 0. The idea is to adhere to the direction of
    %the optimised flux and use the optimised and constrained value as the limit in that
    %direction.
    
    %By default, the flux through the lactase-mediated pathway is stronger
    %than that through the mitochondrial pathaway. Therefore, we will knock
    %down the lactase-mediated pathway and let the mitochondrial pathway
    %respond freely.
    dummy=BN_glycolysis;
    for i = 1:size(Treg_glycolysis)
        if model1.ub(Treg_glycolysis(i))*model1.lb(Treg_glycolysis(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glycolysis(i)), abs(solution.x(Treg_glycolysis(i))*dummy), 'u');
        elseif solution.x(Treg_glycolysis(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glycolysis(i)), solution.x(Treg_glycolysis(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glycolysis(i)), 0, 'l');
        elseif solution.x(Treg_glycolysis(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glycolysis(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glycolysis(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_glycolysis(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glycolysis(i)), solution.x(Treg_glycolysis(i))*dummy, 'l');
        end
    end
    
    dummy=BN_glu_uptake;
    for i = 1:size(Treg_glu_tran)
        if model1.ub(Treg_glu_tran(i))*model1.lb(Treg_glu_tran(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_tran(i)), abs(solution.x(Treg_glu_tran(i))*dummy), 'u');
        elseif solution.x(Treg_glu_tran(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_tran(i)), solution.x(Treg_glu_tran(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_tran(i)), 0, 'l');
        elseif solution.x(Treg_glu_tran(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_tran(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_tran(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_tran(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_tran(i)), solution.x(Treg_glu_tran(i))*dummy, 'l');            
        end
    end
    for i = 1:size(Treg_glu_ex)
        if model1.ub(Treg_glu_ex(i))*model1.lb(Treg_glu_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_ex(i)), abs(solution.x(Treg_glu_ex(i))*dummy), 'u');
        elseif solution.x(Treg_glu_ex(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_ex(i)), solution.x(Treg_glu_ex(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_ex(i)), 0, 'l');
        elseif solution.x(Treg_glu_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_ex(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glu_ex(i)), solution.x(Treg_glu_ex(i))*dummy, 'l');            
        end
    end

    dummy=BN_aa_tran;
    for i = 1:size(Treg_aa_tran)
        if model1.ub(Treg_aa_tran(i))*model1.lb(Treg_aa_tran(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_tran(i)), abs(solution.x(Treg_aa_tran(i))*dummy), 'u');
        elseif solution.x(Treg_aa_tran(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_tran(i)), solution.x(Treg_aa_tran(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_tran(i)), 0, 'l');
        elseif solution.x(Treg_aa_tran(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_tran(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_tran(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_tran(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_tran(i)), solution.x(Treg_aa_tran(i))*dummy, 'l');
            
        end
    end
    for i = 1:size(Treg_aa_ex)
        if model1.ub(Treg_aa_ex(i))*model1.lb(Treg_aa_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_ex(i)), abs(solution.x(Treg_aa_ex(i))*dummy), 'u');
        elseif solution.x(Treg_aa_ex(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_ex(i)), solution.x(Treg_aa_ex(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_ex(i)), 0, 'l');
        elseif solution.x(Treg_aa_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_ex(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_aa_ex(i)), solution.x(Treg_aa_ex(i))*dummy, 'l');
        end
    end
    
    dummy=BN_lip_eff;
    for i = 1:size(Treg_lip_ex)
        if model1.ub(Treg_lip_ex(i))*model1.lb(Treg_lip_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ex(i)), abs(solution.x(Treg_lip_ex(i))*dummy), 'u');
        elseif solution.x(Treg_lip_ex(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ex(i)), solution.x(Treg_lip_ex(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ex(i)), 0, 'l');
        elseif solution.x(Treg_lip_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ex(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ex(i)), solution.x(Treg_lip_ex(i))*dummy, 'l');            
        end
    end
    
    dummy=BN_lip_syn;
    for i = 1:size(Treg_lip_syn)
        if model1.ub(Treg_lip_syn(i))*model1.lb(Treg_lip_syn(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_syn(i)), abs(solution.x(Treg_lip_syn(i))*dummy), 'u');
        elseif solution.x(Treg_lip_syn(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_syn(i)), solution.x(Treg_lip_syn(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_syn(i)), 0, 'l');
        elseif solution.x(Treg_lip_syn(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_syn(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_syn(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_syn(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_syn(i)), solution.x(Treg_lip_syn(i))*dummy, 'l');            
        end
    end
    
    dummy=BN_mit_ox;
    for i = 1:size(Treg_ox_phos)
        if model1.ub(Treg_ox_phos(i))*model1.lb(Treg_ox_phos(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_ox_phos(i)), abs(solution.x(Treg_ox_phos(i))*dummy), 'u');
        elseif solution.x(Treg_ox_phos(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_ox_phos(i)), solution.x(Treg_ox_phos(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_ox_phos(i)), 0, 'l');
        elseif solution.x(Treg_ox_phos(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_ox_phos(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_ox_phos(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_ox_phos(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_ox_phos(i)), solution.x(Treg_ox_phos(i))*dummy, 'l');            
        end
    end
    for i = 1:size(Treg_lip_ox)
        if model1.ub(Treg_lip_ox(i))*model1.lb(Treg_lip_ox(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ox(i)), abs(solution.x(Treg_lip_ox(i))*dummy), 'u');
        elseif solution.x(Treg_lip_ox(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ox(i)), solution.x(Treg_lip_ox(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ox(i)), 0, 'l');
        elseif solution.x(Treg_lip_ox(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ox(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ox(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ox(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_lip_ox(i)), solution.x(Treg_lip_ox(i))*dummy, 'l');            
        end
    end    
    
    dummy=BN_glutaminolysis;
    for i = 1:size(Treg_glutaminolysis)
        if model1.ub(Treg_glutaminolysis(i))*model1.lb(Treg_glutaminolysis(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glutaminolysis(i)), abs(solution.x(Treg_glutaminolysis(i))*dummy), 'u');
        elseif solution.x(Treg_glutaminolysis(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glutaminolysis(i)), solution.x(Treg_glutaminolysis(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glutaminolysis(i)), 0, 'l');
        elseif solution.x(Treg_glutaminolysis(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Treg_glutaminolysis(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glutaminolysis(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Treg_glutaminolysis(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Treg_glutaminolysis(i)), solution.x(Treg_glutaminolysis(i))*dummy, 'l');            
        end
    end
    
    %Set the objective function according to what you want to optimise.
    if Phase==1
        %Set the objective function to biomass minus DNA (G1).
        model1=changeObjective(model1,'Biomass_minusDNA');
        model1=changeRxnBounds(model1, model1.rxns(3862), 0, 'l');
        model1=changeRxnBounds(model1, model1.rxns(3862), 1000, 'u');
        model1=changeRxnBounds(model1, model1.rxns(3855), 0, 'b');
    elseif Phase==2
        %Set the objective function to DNA (S).
        model1=changeObjective(model1,'biomass_DNA');
        model1=changeRxnBounds(model1, model1.rxns(3855), 0, 'l');
        model1=changeRxnBounds(model1, model1.rxns(3855), 1000, 'u');
        model1=changeRxnBounds(model1, model1.rxns(3862), 0, 'b');
    end
    
    %Optimise the modified model.
    solution1=optimizeCbModel(model1);
    
    %Calculate the production rate of ATP.
    ATP=0; %mmol per gram of dry weight per hour.
    for i = 1:size(Treg_ATP_fluxes)
        ATP=ATP+solution1.x(Treg_ATP_fluxes(i))*Treg_ATP_coeff(i);
    end
    
    %Calculate the production rate of AMP.
    AMP=0; %mmol per gram of dry weight per hour.
    for i = 1:size(Treg_AMP_fluxes)
        AMP=AMP+solution1.x(Treg_AMP_fluxes(i))*Treg_AMP_coeff(i);
    end
    
    %Extract the other production rates.
    Biomass=solution1.x(3855); %Growth rate (per hour).
    Biomass_minus_DNA=solution1.x(3862); %Growth rate (per hour).
    DNA=solution1.x(3857); %mmol per gram of dry weight per hour.
    Protein=solution1.x(3856); %mmol per gram of dry weight per hour.
    
    %Combine the production rates.
    output=[Biomass, Biomass_minus_DNA, DNA, Protein, ATP, AMP];
end