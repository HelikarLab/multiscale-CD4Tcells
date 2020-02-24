function output = Th17_MN_control(Th17_model, Phase, solution, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, BN_glycolysis, BN_glu_uptake, BN_aa_tran, BN_mit_ox, BN_lip_eff, BN_glutaminolysis, BN_lip_syn)
    %This function uses the outputs of the BN model to parametrise the MN
    %model for phenotype 14 (Th17) and then calculates various production
    %rates.
    
    %Make a copy of the metabolic model for Th17.
    model1=Th17_model;
    
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
    for i = 1:size(Th17_glycolysis)
        if model1.ub(Th17_glycolysis(i))*model1.lb(Th17_glycolysis(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glycolysis(i)), abs(solution.x(Th17_glycolysis(i))*dummy), 'u');
        elseif solution.x(Th17_glycolysis(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glycolysis(i)), solution.x(Th17_glycolysis(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glycolysis(i)), 0, 'l');
        elseif solution.x(Th17_glycolysis(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glycolysis(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glycolysis(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_glycolysis(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glycolysis(i)), solution.x(Th17_glycolysis(i))*dummy, 'l');
        end
    end
    
    dummy=BN_glu_uptake;
    for i = 1:size(Th17_glu_tran)
        if model1.ub(Th17_glu_tran(i))*model1.lb(Th17_glu_tran(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_tran(i)), abs(solution.x(Th17_glu_tran(i))*dummy), 'u');
        elseif solution.x(Th17_glu_tran(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_tran(i)), solution.x(Th17_glu_tran(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_tran(i)), 0, 'l');
        elseif solution.x(Th17_glu_tran(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_tran(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_tran(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_tran(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_tran(i)), solution.x(Th17_glu_tran(i))*dummy, 'l');            
        end
    end
    for i = 1:size(Th17_glu_ex)
        if model1.ub(Th17_glu_ex(i))*model1.lb(Th17_glu_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_ex(i)), abs(solution.x(Th17_glu_ex(i))*dummy), 'u');
        elseif solution.x(Th17_glu_ex(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_ex(i)), solution.x(Th17_glu_ex(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_ex(i)), 0, 'l');
        elseif solution.x(Th17_glu_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_ex(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glu_ex(i)), solution.x(Th17_glu_ex(i))*dummy, 'l');            
        end
    end

    dummy=BN_aa_tran;
    for i = 1:size(Th17_aa_tran)
        if model1.ub(Th17_aa_tran(i))*model1.lb(Th17_aa_tran(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_tran(i)), abs(solution.x(Th17_aa_tran(i))*dummy), 'u');
        elseif solution.x(Th17_aa_tran(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_tran(i)), solution.x(Th17_aa_tran(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_tran(i)), 0, 'l');
        elseif solution.x(Th17_aa_tran(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_tran(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_tran(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_tran(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_tran(i)), solution.x(Th17_aa_tran(i))*dummy, 'l');
            
        end
    end
    for i = 1:size(Th17_aa_ex)
        if model1.ub(Th17_aa_ex(i))*model1.lb(Th17_aa_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_ex(i)), abs(solution.x(Th17_aa_ex(i))*dummy), 'u');
        elseif solution.x(Th17_aa_ex(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_ex(i)), solution.x(Th17_aa_ex(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_ex(i)), 0, 'l');
        elseif solution.x(Th17_aa_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_ex(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_aa_ex(i)), solution.x(Th17_aa_ex(i))*dummy, 'l');
        end
    end
    
    dummy=BN_lip_eff;
    for i = 1:size(Th17_lip_ex)
        if model1.ub(Th17_lip_ex(i))*model1.lb(Th17_lip_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ex(i)), abs(solution.x(Th17_lip_ex(i))*dummy), 'u');
        elseif solution.x(Th17_lip_ex(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ex(i)), solution.x(Th17_lip_ex(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ex(i)), 0, 'l');
        elseif solution.x(Th17_lip_ex(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ex(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ex(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ex(i)), solution.x(Th17_lip_ex(i))*dummy, 'l');            
        end
    end
    
    dummy=BN_lip_syn;
    for i = 1:size(Th17_lip_syn)
        if model1.ub(Th17_lip_syn(i))*model1.lb(Th17_lip_syn(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_syn(i)), abs(solution.x(Th17_lip_syn(i))*dummy), 'u');
        elseif solution.x(Th17_lip_syn(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_syn(i)), solution.x(Th17_lip_syn(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_syn(i)), 0, 'l');
        elseif solution.x(Th17_lip_syn(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_syn(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_syn(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_syn(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_syn(i)), solution.x(Th17_lip_syn(i))*dummy, 'l');            
        end
    end
    
    dummy=BN_mit_ox;
    for i = 1:size(Th17_ox_phos)
        if model1.ub(Th17_ox_phos(i))*model1.lb(Th17_ox_phos(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_ox_phos(i)), abs(solution.x(Th17_ox_phos(i))*dummy), 'u');
        elseif solution.x(Th17_ox_phos(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_ox_phos(i)), solution.x(Th17_ox_phos(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_ox_phos(i)), 0, 'l');
        elseif solution.x(Th17_ox_phos(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_ox_phos(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_ox_phos(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_ox_phos(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_ox_phos(i)), solution.x(Th17_ox_phos(i))*dummy, 'l');            
        end
    end
    for i = 1:size(Th17_lip_ox)
        if model1.ub(Th17_lip_ox(i))*model1.lb(Th17_lip_ox(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ox(i)), abs(solution.x(Th17_lip_ox(i))*dummy), 'u');
        elseif solution.x(Th17_lip_ox(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ox(i)), solution.x(Th17_lip_ox(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ox(i)), 0, 'l');
        elseif solution.x(Th17_lip_ox(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ox(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ox(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ox(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_lip_ox(i)), solution.x(Th17_lip_ox(i))*dummy, 'l');            
        end
    end    
    
    dummy=BN_glutaminolysis;
    for i = 1:size(Th17_glutaminolysis)
        if model1.ub(Th17_glutaminolysis(i))*model1.lb(Th17_glutaminolysis(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glutaminolysis(i)), abs(solution.x(Th17_glutaminolysis(i))*dummy), 'u');
        elseif solution.x(Th17_glutaminolysis(i))>0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glutaminolysis(i)), solution.x(Th17_glutaminolysis(i))*dummy, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glutaminolysis(i)), 0, 'l');
        elseif solution.x(Th17_glutaminolysis(i))==0
            model1=changeRxnBounds(model1, model1.rxns(Th17_glutaminolysis(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glutaminolysis(i)), 0, 'l');
        else
            model1=changeRxnBounds(model1, model1.rxns(Th17_glutaminolysis(i)), 0, 'u');
            model1=changeRxnBounds(model1, model1.rxns(Th17_glutaminolysis(i)), solution.x(Th17_glutaminolysis(i))*dummy, 'l');            
        end
    end
    
    %Set the objective function according to what you want to optimise.
    if Phase==1
        %Set the objective function to biomass minus DNA (G1).
        model1=changeObjective(model1,'Biomass_minusDNA');
        model1=changeRxnBounds(model1, model1.rxns(5231), 0, 'l');
        model1=changeRxnBounds(model1, model1.rxns(5231), 1000, 'u');
        model1=changeRxnBounds(model1, model1.rxns(5224), 0, 'b');
    elseif Phase==2
        %Set the objective function to DNA (S).
        model1=changeObjective(model1,'biomass_DNA');
        model1=changeRxnBounds(model1, model1.rxns(5224), 0, 'l');
        model1=changeRxnBounds(model1, model1.rxns(5224), 1000, 'u');
        model1=changeRxnBounds(model1, model1.rxns(5231), 0, 'b');
    end
    
    %Optimise the modified model.
    solution1=optimizeCbModel(model1);
    
    %Calculate the production rate of ATP.
    ATP=0; %mmol per gram of dry weight per hour.
    for i = 1:size(Th17_ATP_fluxes)
        ATP=ATP+solution1.x(Th17_ATP_fluxes(i))*Th17_ATP_coeff(i);
    end
    
    %Calculate the production rate of AMP.
    AMP=0; %mmol per gram of dry weight per hour.
    for i = 1:size(Th17_AMP_fluxes)
        AMP=AMP+solution1.x(Th17_AMP_fluxes(i))*Th17_AMP_coeff(i);
    end
    
    %Extract the other production rates.
    Biomass=solution1.x(5224); %Growth rate (per hour).
    Biomass_minus_DNA=solution1.x(5231); %Growth rate (per hour).
    DNA=solution1.x(5226); %mmol per gram of dry weight per hour.
    Protein=solution1.x(5225); %mmol per gram of dry weight per hour.
    
    %Combine the production rates.
    output=[Biomass, Biomass_minus_DNA, DNA, Protein, ATP, AMP];
end