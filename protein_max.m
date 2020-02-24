function protein_max = protein_max(phenotype, model)
    %This function calculates the maximum protein synthesis rate for each
    %major phenotype.

    %Calculation: maximum protein synthesis rate (phenotype 0, Th0).
    if phenotype == 0
        model=changeObjective(model,'biomass_protein');
        model=changeRxnBounds(model, model.rxns(4235), 0, 'l');
        model=changeRxnBounds(model, model.rxns(4235), 1000, 'u');
        model=changeRxnBounds(model, model.rxns(4242), 0, 'b');
        solution=optimizeCbModel(model);
        protein_max=solution.x(4236); %mmol per gram of dry weight per hour.
    end
    
    %Calculation: maximum protein synthesis rate (phenotype 12, Th1).
    if phenotype == 12
        model=changeObjective(model,'biomass_protein');
        model=changeRxnBounds(model, model.rxns(4161), 0, 'l');
        model=changeRxnBounds(model, model.rxns(4161), 1000, 'u');
        model=changeRxnBounds(model, model.rxns(4168), 0, 'b');
        solution=optimizeCbModel(model);
        protein_max=solution.x(4162); %mmol per gram of dry weight per hour.
    end
    
    %Calculation: maximum protein synthesis rate (phenotype 13, Th2).
    if phenotype == 13
        model=changeObjective(model,'biomass_protein');
        model=changeRxnBounds(model, model.rxns(4675), 0, 'l');
        model=changeRxnBounds(model, model.rxns(4675), 1000, 'u');
        model=changeRxnBounds(model, model.rxns(4682), 0, 'b');
        solution=optimizeCbModel(model);
        protein_max=solution.x(4676); %mmol per gram of dry weight per hour.
    end
    
    %Calculation: maximum protein synthesis rate (phenotype 14, Th17).
    if phenotype == 14
        model=changeObjective(model,'biomass_protein');
        model=changeRxnBounds(model, model.rxns(5224), 0, 'l');
        model=changeRxnBounds(model, model.rxns(5224), 1000, 'u');
        model=changeRxnBounds(model, model.rxns(5231), 0, 'b');
        solution=optimizeCbModel(model);
        protein_max=solution.x(5225); %mmol per gram of dry weight per hour.
    end
    
    %Calculation: maximum protein synthesis rate (phenotype 15, Treg).
    if phenotype == 15
        model=changeObjective(model,'biomass_protein');
        model=changeRxnBounds(model, model.rxns(3855), 0, 'l');
        model=changeRxnBounds(model, model.rxns(3855), 1000, 'u');
        model=changeRxnBounds(model, model.rxns(3862), 0, 'b');
        solution=optimizeCbModel(model);
        protein_max=solution.x(3856); %mmol per gram of dry weight per hour.
    end

end