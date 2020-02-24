%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Timekeeping starts.
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialise supporting packages and load supplementary files.

%Load Cobra toolbox.
initCobraToolbox;

%Population dynamics of the antigen population.
load('antigen.mat');

%Files for Th0 metabolism (phenotype 0).
load('Th0_model.mat');
load('Th0_table.mat');
load('Th0_MNlibrary.mat');

%Files for Th1 metabolism (phenotype 12).
load('Th1_model_v1.8.mat');
load('Th1_table.mat');
load('Th1_MNlibrary.mat');

%Files for Th2 metabolism (phenotype 13).
load('Th2_model.mat');
load('Th2_table.mat');
load('Th2_MNlibrary.mat');

%Files for Th17 metabolism (phenotype 14).
load('Th17_model.mat');
load('Th17_table.mat');
load('Th17_MNlibrary.mat');

%Files for Treg metabolism (phenotype 15).
load('Treg_Model.mat');
load('Treg_table.mat');
load('Treg_MNlibrary.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Decide on the duration of the simulation by indicating the
%number of time steps and the size of each step (hours). Also, decide how
%many times the model should be simulated.

%According to a paper, titled 'Regulation of CD4+ T-cell Contraction During
%Pathogen Challenge' and authored by McKinstry in 2010, the viral load of
%influenza peaks 2 days after infection and begins to fall 3 days later,
%taking 5 days to clear.
steps=400;
step_size=1;

%This should be fixed by a convergence study.
rounds=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Provide the population dynamics of the antigen under consideration.

%According to a paper, titled 'Regulation of CD4+ T-cell Contraction During
%Pathogen Challenge' and authored by McKinstry in 2010, the viral load of
%influenza peaks 2 days after infection and begins to fall 3 days later,
%taking 5 days to clear.
TCR=zeros(3, steps+1);
TCR(1, :)=Antigen_target;
TCR(2, :)=Antigen_LN;
TCR(3, :)=Antigen_circulation;
CD28=zeros(3, steps+1);
CD28(1, :)=Antigen_target;
CD28(2, :)=Antigen_LN;
CD28(3, :)=Antigen_circulation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the population of CD4+ T cells.

%Initial conditions: numbers of CD4+ T cells in different compartments (PATHOGEN-DEPENDENT).
%1-100 cells per million for each epitope. CD4 T-cell Memory Generation and Maintenance (Gasper 2014).
%175 influenza A epitopes for CD4+ T cells. Ab and T Cell Epitopes of Influenza A Virus, Knowledge and Opportunities (Bui 2006).
%According to Alyssa's document about Q2, there are 1.69e8*(175*50/1e6)=1.48e6 cells in the lungs, 3.13e10/600*(175*50/1e6)=4.56e5 cells in a lymph node (600 lymph nodes in a human), and 2.04e9*(175*50/1e6)=1.785e7 cells in blood and lymph.
Nt0=32; %Representative cell count in T.
Nln0=10; %Representative cell count in LN.
Nc0=391; %Representative cell count in C.
Nt0_all=1.48e6; %Actual cell count in T.
Nln0_all=4.56e5; %Actual cell count in LN.
Nc0_all=1.785e7; %Actual cell count in C.

k_effector_CtoT=0.02083*step_size; %An effector cell exits LN and travels in C to reach T. Harnessing the Heterogeneity of T Cell Differentiation Fate to Fine-tune Generation of Effector and Memory T Cells (Gong 2014). Regulation of CD4+ T-cell Contraction During Pathogen Challenge (McKinstry 2010).
act_thres_mem=2/step_size; %Number of TCR stimuli required to reactivate a memory T cell. T Cell Receptor and Cytokine Signaling Can Function at Different Stages to Establish and Maintain Transcriptional Memory and Enable T Helper Cell Differentiation (Bevington 2017).
k_toln_boost=5; %When TCR signalling is on in LN, and if the naive cell is in T, the probability of cell migration to LN increases by five times. Innate Control of Adaptive Immunity via Remodeling of Lymph Node Feed Arteriole (Soderberg 2005).
TCR_max=24/step_size; %The maximum TCR signalling strength. Two Distinct Stages in the Transition from Naive CD4 T Cells to Effectors, Early Antigen-dependent and Late Cytokine-driven Expansion and Differentiation (Jelley-Gibbs 2000).
CD28_max=24/step_size; %The maximum CD28 signalling strength. Two Distinct Stages in the Transition from Naive CD4 T Cells to Effectors, Early Antigen-dependent and Late Cytokine-driven Expansion and Differentiation (Jelley-Gibbs 2000).
Ease_restim=0.5; %This parameter increases with the ease of restimulating an effector cell or reactivated memory cell. Very useful for calibration.
Ease_ACAD=0.5; %This parameter increases with the ease with which an effector cell or reactivated memory cell dies due to ACAD. It is useful for calibration.
Ease_memory=0.05; %This parameter increases with the ease with which an effector cell or reactivated memory cell forms memory. It is useful for calibration.
Ease_relax=0.5; %This parameter increases with the ease with which an effector cell or reactivated memory cell relaxes in the absence of restimulation. It is useful for calibration.
Refill_thres=0.2; %When the naive cell population drops below this percentage of the initial population, naive cells are added to the system.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialise and parametrise the dynamics of cytokines.

%Parameters: compartment volumes and flow rates (PATHOGEN-DEPENDENT).
Qa=356; %Blood flow rate into the lungs (L per hour). A Bottom-Up Whole-Body Physiologically Based Pharmacokinetic Model to Mechanistically Predict Tissue Distribution and the Rate of Subcutaneous Absorption of Therapeutic Proteins (Gill 2016).
Qb=0.312*0.03; %Lymph flow rate into the draining lymph node (L per hour). A Bottom-Up Whole-Body Physiologically Based Pharmacokinetic Model to Mechanistically Predict Tissue Distribution and the Rate of Subcutaneous Absorption of Therapeutic Proteins (Gill 2016).
Vt=0.547; %Volume of the target organ, the lungs (L). A Bottom-Up Whole-Body Physiologically Based Pharmacokinetic Model to Mechanistically Predict Tissue Distribution and the Rate of Subcutaneous Absorption of Therapeutic Proteins (Gill 2016).
Vln=4.545e-4; %Volume of the draining lymph node (L). Normal mediastinal lymph node size and number: CT and anatomic study (Genereux 1984). Size and volume of metastatic and non-metastatic lymph nodes in pelvis and lower abdomen in patients with carcinoma of the bladder undergoing radical cystectomy (Jensen 2010).
Vc=3.5+0.3; %Volume of the circulation system: blood plus lymph (L). A Bottom-Up Whole-Body Physiologically Based Pharmacokinetic Model to Mechanistically Predict Tissue Distribution and the Rate of Subcutaneous Absorption of Therapeutic Proteins (Gill 2016).

%Parameters: cytokine production rates due to the pathogen (PATHOGEN-DEPENDENT).
PIL2in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIL4in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIL6in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIL12in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIL17in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIL18in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIL21in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIL23in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIL27in=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PIFNgin=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.
PTGFbin=[0, 0, 0]; %Pathogen-induced production rate (M per hour). Excel file.

%Parameters and initial conditions: cytokine dynamics (FIXED).
IL20=[2.83639664170638E-10, 2.83639664170638E-10, 2.83639664170638E-10]; %Initial concentrations (M). Excel file.
IL40=[1.14337983077979E-09, 1.14337983077979E-09, 1.14337983077979E-09]; %Initial concentrations (M). Excel file.
IL60=[8.43241420018551E-10, 8.43241420018551E-10, 8.43241420018551E-10]; %Initial concentrations (M). Excel file.
IL120=[1.61256510731621E-10, 1.61256510731621E-10, 1.61256510731621E-10]; %Initial concentrations (M). Excel file.
IL170=[0, 0, 0]; %Initial concentrations (M). Excel file.
IL180=[0, 0, 0]; %Initial concentrations (M). Excel file.
IL210=[0, 0, 0]; %Initial concentrations (M). Excel file.
IL230=[0, 0, 0]; %Initial concentrations (M). Excel file.
IL270=[0, 0, 0]; %Initial concentrations (M). Excel file.
IFNg0=[2.58424643373992E-10, 2.58424643373992E-10, 2.58424643373992E-10]; %Initial concentrations (M). Excel file.
TGFb0=[1.12762454613112E-10, 1.12762454613112E-10, 1.12762454613112E-10]; %Initial concentrations (M). Excel file.
tau=step_size; %Use the step size (hours) as the time scale.
kdeg=0.1; %Degradation rate constant (per hour). Mathematical Modelling of Immune Response in Tissues (Su 2009).
Cs=[1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8]; %Concentration scales (M). Cytokines as Biomarkers in Rheumatoid Arthritis (Burska 2014).

%Calculation: cytokine production rates required to make the initial cytokine concentrations the baseline concentrations.
PIL2=[IL20(1)*kdeg-Qa*IL20(3)/Vt+Qb*IL20(1)/Vt+(Qa-Qb)*IL20(1)/Vt, IL20(2)*kdeg-Qb*IL20(1)/Vln+Qb*IL20(2)/Vln, kdeg*IL20(3)-(Qa-Qb)*IL20(1)/Vc-Qb*IL20(2)/Vc+Qa*IL20(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIL4=[IL40(1)*kdeg-Qa*IL40(3)/Vt+Qb*IL40(1)/Vt+(Qa-Qb)*IL40(1)/Vt, IL40(2)*kdeg-Qb*IL40(1)/Vln+Qb*IL40(2)/Vln, kdeg*IL40(3)-(Qa-Qb)*IL40(1)/Vc-Qb*IL40(2)/Vc+Qa*IL40(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIL6=[IL60(1)*kdeg-Qa*IL60(3)/Vt+Qb*IL60(1)/Vt+(Qa-Qb)*IL60(1)/Vt, IL60(2)*kdeg-Qb*IL60(1)/Vln+Qb*IL60(2)/Vln, kdeg*IL60(3)-(Qa-Qb)*IL60(1)/Vc-Qb*IL60(2)/Vc+Qa*IL60(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIL12=[IL120(1)*kdeg-Qa*IL120(3)/Vt+Qb*IL120(1)/Vt+(Qa-Qb)*IL120(1)/Vt, IL120(2)*kdeg-Qb*IL120(1)/Vln+Qb*IL120(2)/Vln, kdeg*IL120(3)-(Qa-Qb)*IL120(1)/Vc-Qb*IL120(2)/Vc+Qa*IL120(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIL17=[IL170(1)*kdeg-Qa*IL170(3)/Vt+Qb*IL170(1)/Vt+(Qa-Qb)*IL170(1)/Vt, IL170(2)*kdeg-Qb*IL170(1)/Vln+Qb*IL170(2)/Vln, kdeg*IL170(3)-(Qa-Qb)*IL170(1)/Vc-Qb*IL170(2)/Vc+Qa*IL170(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIL18=[IL180(1)*kdeg-Qa*IL180(3)/Vt+Qb*IL180(1)/Vt+(Qa-Qb)*IL180(1)/Vt, IL180(2)*kdeg-Qb*IL180(1)/Vln+Qb*IL180(2)/Vln, kdeg*IL180(3)-(Qa-Qb)*IL180(1)/Vc-Qb*IL180(2)/Vc+Qa*IL180(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIL21=[IL210(1)*kdeg-Qa*IL210(3)/Vt+Qb*IL210(1)/Vt+(Qa-Qb)*IL210(1)/Vt, IL210(2)*kdeg-Qb*IL210(1)/Vln+Qb*IL210(2)/Vln, kdeg*IL210(3)-(Qa-Qb)*IL210(1)/Vc-Qb*IL210(2)/Vc+Qa*IL210(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIL23=[IL230(1)*kdeg-Qa*IL230(3)/Vt+Qb*IL230(1)/Vt+(Qa-Qb)*IL230(1)/Vt, IL230(2)*kdeg-Qb*IL230(1)/Vln+Qb*IL230(2)/Vln, kdeg*IL230(3)-(Qa-Qb)*IL230(1)/Vc-Qb*IL230(2)/Vc+Qa*IL230(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIL27=[IL270(1)*kdeg-Qa*IL270(3)/Vt+Qb*IL270(1)/Vt+(Qa-Qb)*IL270(1)/Vt, IL270(2)*kdeg-Qb*IL270(1)/Vln+Qb*IL270(2)/Vln, kdeg*IL270(3)-(Qa-Qb)*IL270(1)/Vc-Qb*IL270(2)/Vc+Qa*IL270(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PIFNg=[IFNg0(1)*kdeg-Qa*IFNg0(3)/Vt+Qb*IFNg0(1)/Vt+(Qa-Qb)*IFNg0(1)/Vt, IFNg0(2)*kdeg-Qb*IFNg0(1)/Vln+Qb*IFNg0(2)/Vln, kdeg*IFNg0(3)-(Qa-Qb)*IFNg0(1)/Vc-Qb*IFNg0(2)/Vc+Qa*IFNg0(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.
PTGFb=[TGFb0(1)*kdeg-Qa*TGFb0(3)/Vt+Qb*TGFb0(1)/Vt+(Qa-Qb)*TGFb0(1)/Vt, TGFb0(2)*kdeg-Qb*TGFb0(1)/Vln+Qb*TGFb0(2)/Vln, kdeg*TGFb0(3)-(Qa-Qb)*TGFb0(1)/Vc-Qb*TGFb0(2)/Vc+Qa*TGFb0(3)/Vc]; %Production rates (M per hour) in different compartments (T, LN, C). Excel file.

%Calculation: cytokine production rates in CD4+ T cells.
PIL2ex=[10*3600*Nt0_all/Nt0/6.02e23/Vt, 10*3600*Nln0_all/Nln0/6.02e23/Vln, 10*3600*Nc0_all/Nc0/6.02e23/Vc]; %Production rate due to CD4+ T cells (M per hour). Polyfunctional responses by human T cells result from sequential release of cytokines (Han 2011). Excel file.
PIL4ex=[10*3600*Nt0_all/Nt0/6.02e23/Vt, 10*3600*Nln0_all/Nln0/6.02e23/Vln, 10*3600*Nc0_all/Nc0/6.02e23/Vc]; %Production rate due to CD4+ T cells (M per hour). Polyfunctional responses by human T cells result from sequential release of cytokines (Han 2011). Excel file.
PIL6ex=[10*3600*Nt0_all/Nt0/6.02e23/Vt, 10*3600*Nln0_all/Nln0/6.02e23/Vln, 10*3600*Nc0_all/Nc0/6.02e23/Vc]; %Production rate due to CD4+ T cells (M per hour). Polyfunctional responses by human T cells result from sequential release of cytokines (Han 2011). Excel file.
PIL12ex=[0, 0, 0]; %Production rate due to CD4+ T cells (M per hour). Excel file.
PIL17ex=[10*3600*Nt0_all/Nt0/6.02e23/Vt, 10*3600*Nln0_all/Nln0/6.02e23/Vln, 10*3600*Nc0_all/Nc0/6.02e23/Vc]; %Production rate due to CD4+ T cells (M per hour). Polyfunctional responses by human T cells result from sequential release of cytokines (Han 2011). Excel file.
PIL18ex=[0, 0, 0]; %Production rate due to CD4+ T cells (M per hour). Excel file.
PIL21ex=[10*3600*Nt0_all/Nt0/6.02e23/Vt, 10*3600*Nln0_all/Nln0/6.02e23/Vln, 10*3600*Nc0_all/Nc0/6.02e23/Vc]; %Production rate due to CD4+ T cells (M per hour). Polyfunctional responses by human T cells result from sequential release of cytokines (Han 2011). Excel file.
PIL23ex=[0, 0, 0]; %Production rate due to CD4+ T cells (M per hour). Excel file.
PIL27ex=[0, 0, 0]; %Production rate due to CD4+ T cells (M per hour). Excel file.
PIFNgex=[10*3600*Nt0_all/Nt0/6.02e23/Vt, 10*3600*Nln0_all/Nln0/6.02e23/Vln, 10*3600*Nc0_all/Nc0/6.02e23/Vc]; %Production rate due to CD4+ T cells (M per hour). Polyfunctional responses by human T cells result from sequential release of cytokines (Han 2011). Excel file.
PTGFbex=[0, 0, 0]; %Production rate due to CD4+ T cells (M per hour). Excel file.

%Calculation: nondimensionalisation of the initial cytokine concentrations.
IL20=IL20/Cs(1); %Nondimensionalised initial concentrations.
IL40=IL40/Cs(2); %Nondimensionalised initial concentrations.
IL60=IL60/Cs(3); %Nondimensionalised initial concentrations.
IL120=IL120/Cs(4); %Nondimensionalised initial concentrations.
IL170=IL170/Cs(5); %Nondimensionalised initial concentrations.
IL180=IL180/Cs(6); %Nondimensionalised initial concentrations.
IL210=IL210/Cs(7); %Nondimensionalised initial concentrations.
IL230=IL230/Cs(8); %Nondimensionalised initial concentrations.
IL270=IL270/Cs(9); %Nondimensionalised initial concentrations.
IFNg0=IFNg0/Cs(10); %Nondimensionalised initial concentrations.
TGFb0=TGFb0/Cs(11); %Nondimensionalised initial concentrations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the default metabolic fluxes for the five major phenotypes.

%Calculate the maximum protein synthesis rate for each major phenotype.
protein_max_vector=zeros(1, 5); %mmol per gram of dry weight per hour.
protein_max_vector(1, 1)=protein_max(0, Th0); %Phenotype 0 (Th0).
protein_max_vector(1, 2)=protein_max(12, Th1); %Phenotype 12 (Th1).
protein_max_vector(1, 3)=protein_max(13, Th2); %Phenotype 13 (Th2).
protein_max_vector(1, 4)=protein_max(14, Th17); %Phenotype 14 (Th17).
protein_max_vector(1, 5)=protein_max(15, Treg); %Phenotype 15 (Treg).

%Optimise the metabolic model for phenotype 0 (Th0). Need to do it twice, once for G1 and once for S.
model0=Th0;
%Set the objective function to biomass minus DNA (G1).
model0=changeObjective(model0,'Biomass_minusDNA');
model0=changeRxnBounds(model0, model0.rxns(4242), 0, 'l');
model0=changeRxnBounds(model0, model0.rxns(4242), 1000, 'u');
model0=changeRxnBounds(model0, model0.rxns(4235), 0, 'b');
solution0_1=optimizeCbModel(model0);
%Set the objective function to DNA (S).
model0=changeObjective(model0,'biomass_DNA');
model0=changeRxnBounds(model0, model0.rxns(4235), 0, 'l');
model0=changeRxnBounds(model0, model0.rxns(4235), 1000, 'u');
model0=changeRxnBounds(model0, model0.rxns(4242), 0, 'b');
solution0_2=optimizeCbModel(model0);

%Optimise the metabolic model for phenotype 12 (Th1). Need to do it twice, once for G1 and once for S.
model12=Th1;
%Set the objective function to biomass minus DNA (G1).
model12=changeObjective(model12,'Biomass_minusDNA');
model12=changeRxnBounds(model12, model12.rxns(4168), 0, 'l');
model12=changeRxnBounds(model12, model12.rxns(4168), 1000, 'u');
model12=changeRxnBounds(model12, model12.rxns(4161), 0, 'b');
solution12_1=optimizeCbModel(model12);
%Set the objective function to DNA (S).
model12=changeObjective(model12,'biomass_DNA');
model12=changeRxnBounds(model12, model12.rxns(4161), 0, 'l');
model12=changeRxnBounds(model12, model12.rxns(4161), 1000, 'u');
model12=changeRxnBounds(model12, model12.rxns(4168), 0, 'b');
solution12_2=optimizeCbModel(model12);

%Optimise the metabolic model for phenotype 13 (Th2). Need to do it twice, once for G1 and once for S.
model13=Th2;
%Set the objective function to biomass minus DNA (G1).
model13=changeObjective(model13,'Biomass_minusDNA');
model13=changeRxnBounds(model13, model13.rxns(4682), 0, 'l');
model13=changeRxnBounds(model13, model13.rxns(4682), 1000, 'u');
model13=changeRxnBounds(model13, model13.rxns(4675), 0, 'b');
solution13_1=optimizeCbModel(model13);
%Set the objective function to DNA (S).
model13=changeObjective(model13,'biomass_DNA');
model13=changeRxnBounds(model13, model13.rxns(4675), 0, 'l');
model13=changeRxnBounds(model13, model13.rxns(4675), 1000, 'u');
model13=changeRxnBounds(model13, model13.rxns(4682), 0, 'b');
solution13_2=optimizeCbModel(model13);

%Optimise the metabolic model for phenotype 14 (Th17). Need to do it twice, once for G1 and once for S.
model14=Th17;
%Set the objective function to biomass minus DNA (G1).
model14=changeObjective(model14,'Biomass_minusDNA');
model14=changeRxnBounds(model14, model14.rxns(5231), 0, 'l');
model14=changeRxnBounds(model14, model14.rxns(5231), 1000, 'u');
model14=changeRxnBounds(model14, model14.rxns(5224), 0, 'b');
solution14_1=optimizeCbModel(model14);
%Set the objective function to DNA (S).
model14=changeObjective(model14,'biomass_DNA');
model14=changeRxnBounds(model14, model14.rxns(5224), 0, 'l');
model14=changeRxnBounds(model14, model14.rxns(5224), 1000, 'u');
model14=changeRxnBounds(model14, model14.rxns(5231), 0, 'b');
solution14_2=optimizeCbModel(model14);

%Optimise the metabolic model for phenotype 15 (Treg). Need to do it twice, once for G1 and once for S.
model15=Treg;
%Set the objective function to biomass minus DNA (G1).
model15=changeObjective(model15,'Biomass_minusDNA');
model15=changeRxnBounds(model15, model15.rxns(3862), 0, 'l');
model15=changeRxnBounds(model15, model15.rxns(3862), 1000, 'u');
model15=changeRxnBounds(model15, model15.rxns(3855), 0, 'b');
solution15_1=optimizeCbModel(model15);
%Set the objective function to DNA (S).
model15=changeObjective(model15,'biomass_DNA');
model15=changeRxnBounds(model15, model15.rxns(3855), 0, 'l');
model15=changeRxnBounds(model15, model15.rxns(3855), 1000, 'u');
model15=changeRxnBounds(model15, model15.rxns(3862), 0, 'b');
solution15_2=optimizeCbModel(model15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run multiple simulations due to the inherently stochastic nature of the model.

%For each simulation, do the following.
for name=1:rounds
    
    %Display the current simulation number.
    name

    %Create arrays to hold the initial cytokine concentrations.
    IL2=IL20;
    IL4=IL40;
    IL6=IL60;
    IL12=IL120;
    IL17=IL170;
    IL18=IL180;
    IL21=IL210;
    IL23=IL230;
    IL27=IL270;
    IFNg=IFNg0;
    TGFb=TGFb0;

    %Create a cell array to hold three populations of naive CD4+ T cells, one for each compartment.
    clear now;
    now=cell(Nt0+Nln0+Nc0, 1);
    
    %Start a new cell ID system for the current simulation.
    ID=0;

    %Create naive CD4+ T cells to populate the cell array.
    for c=1:Nt0
        ID=ID+1;
        now{c}=CD4(ID, 1, step_size, 1-Nln0/Nt0, Nln0/Nt0/2, Nln0/Nt0/2);
    end
    for c=Nt0+1:Nt0+Nln0
        ID=ID+1;
        now{c}=CD4(ID, 2, step_size, 0, 0.5, 0.5);
    end
    for c=Nt0+Nln0+1:Nt0+Nln0+Nc0
        ID=ID+1;
        now{c}=CD4(ID, 3, step_size, Nln0/Nc0, 0, 1-Nln0/Nc0);
    end

    %Create arrays to record the numbers of cytokine-producing CD4+ T cells
    %in the three compartments (T, LN, and C). Initially, the naive CD4+ T
    %cells do not produce any cytokines.
    N_IFNg=[0, 0, 0];
    N_IL4=[0, 0, 0];
    N_IL6=[0, 0, 0];
    N_IL2=[0, 0, 0];
    N_IL17=[0, 0, 0];
    N_IL21=[0, 0, 0];

    %Record the initial conditions: the entire now array and the cytokine concentrations.
    %In log_cytokines, entries 1-11 are the concentrations in T; 12-22, LN; 23-33, C.
    clear log_cytokines;
    cytokines_t=[IL2(1)*Cs(1), IL4(1)*Cs(2), IL6(1)*Cs(3), IL12(1)*Cs(4), IL17(1)*Cs(5), IL18(1)*Cs(6), IL21(1)*Cs(7), IL23(1)*Cs(8), IL27(1)*Cs(9), IFNg(1)*Cs(10), TGFb(1)*Cs(11)];
    cytokines_ln=[IL2(2)*Cs(1), IL4(2)*Cs(2), IL6(2)*Cs(3), IL12(2)*Cs(4), IL17(2)*Cs(5), IL18(2)*Cs(6), IL21(2)*Cs(7), IL23(2)*Cs(8), IL27(2)*Cs(9), IFNg(2)*Cs(10), TGFb(2)*Cs(11)];
    cytokines_c=[IL2(3)*Cs(1), IL4(3)*Cs(2), IL6(3)*Cs(3), IL12(3)*Cs(4), IL17(3)*Cs(5), IL18(3)*Cs(6), IL21(3)*Cs(7), IL23(3)*Cs(8), IL27(3)*Cs(9), IFNg(3)*Cs(10), TGFb(3)*Cs(11)];
    log_cytokines=[cytokines_t, cytokines_ln, cytokines_c];
    fullname=['Simulation', num2str(name), '_0', '.mat'];
    save(fullname, 'now', 'log_cytokines');

    %Temporal dynamics of the current simulation.
    for t=1:steps
        
        %Display the current time step which lasts from time t-1 to time t.
        t

        %Evaluate each cell at the beginning of the current time step.
        %Note that any cells created during this time step will not be
        %evaluated.
        for c=1:size(now, 1)
            
            %If the cell is dead, do nothing.
            if now{c}.State==0
                continue
                
            %If the cell is naive, update its migration probabilities according
            %to the TCR level in LN and stimulate/activate the cell.
            elseif now{c}.State==1
                
                %When TCR signalling is on in LN, and if the naive cell is in
                %T, the probability of cell migration to LN increases by five times.
                %Innate Control of Adaptive Immunity via Remodeling of Lymph Node Feed Arteriole (Soderberg 2005).
                %Assume that the naive cells in T no longer move into C under these circumstances (Assumption).
                if rand<TCR(2, t) & now{c}.Location==1
                    now{c}.k_toln=k_toln_boost*Nln0/Nt0/2;
                    now{c}.k_toc=0;
                    if now{c}.k_toln<1
                        now{c}.k_tot=1-now{c}.k_toln;
                    else
                        now{c}.k_tot=0;
                    end
                %If TCR signalling is on in LN and the naive cell is already in LN, it stays there (assumption).
                elseif rand<TCR(2, t) & now{c}.Location==2
                    now{c}.k_tot=0;
                    now{c}.k_toln=1;
                    now{c}.k_toc=0;
                %Otherwise, it retains or regains the default migration probabilities. 
                else
                    if now{c}.Location==1
                        now{c}.k_tot=1-Nln0/Nt0;
                        now{c}.k_toln=Nln0/Nt0/2;
                        now{c}.k_toc=Nln0/Nt0/2;
                    elseif now{c}.Location==2
                        now{c}.k_tot=0;
                        now{c}.k_toln=0.5;
                        now{c}.k_toc=0.5;
                    else
                        now{c}.k_tot=Nln0/Nc0;
                        now{c}.k_toln=0;
                        now{c}.k_toc=1-Nln0/Nc0;
                    end
                end
                
                %When TCR signalling is on in LN, if the cell is in LN, it is stimulated once and the TCR signalling and CD28 signalling strengths are recorded. The cytokine signals are recorded too.
                %Naive CD4+ T cells activate in LN only. CD4+ T Cell Priming as Biomarker to Study Immune Response to Preventive Vaccines (Ciabattini, 2013).
                if rand<TCR(2, t) & now{c}.Location==2
                    now{c}.TCR_Stimuli=now{c}.TCR_Stimuli+1;
                    now{c}.TCR_Strength=now{c}.TCR_Strength+TCR(2, t);
                    now{c}.CD28_Strength=now{c}.CD28_Strength+CD28(2, t);
                    now{c}.IL2_sum=now{c}.IL2_sum+IL2(2);
                    now{c}.IL4_sum=now{c}.IL4_sum+IL4(2);
                    now{c}.IL6_sum=now{c}.IL6_sum+IL6(2);
                    now{c}.IL12_sum=now{c}.IL12_sum+IL12(2);
                    now{c}.IL18_sum=now{c}.IL18_sum+IL18(2);
                    now{c}.IL23_sum=now{c}.IL23_sum+IL23(2);
                    now{c}.IL27_sum=now{c}.IL27_sum+IL27(2);
                    now{c}.IFNg_sum=now{c}.IFNg_sum+IFNg(2);
                    now{c}.TGFb_sum=now{c}.TGFb_sum+TGFb(2);
                    %Following the said stimulation, if the cell has been stimulated sufficiently, it activates.
                    %Stimulation during activation is sufficient for later expansion and effector functions. Two Distinct Stages in the Transition from Naive CD4 T Cells to Effectors, Early Antigen-dependent and Late Cytokine-driven Expansion and Differentiation (Jelley-Gibbs 2000).
                    %Before running the BN and MN models, we do not know the protein synthesis rate, so we will link it to the TCR strength.
                    %Run the BN model and MN model to decide the activated cell's internal state. 
                    %Finally, update the migration probabilities. Before six divisions, an effector cell cannot exit LN. CD4 Effector T Cell Subsets in the Response to Influenza (Roman 2002).
                    if now{c}.TCR_Stimuli==now{c}.act_thres
                        now{c}.State=2;
                        now{c}.Ribosome=now{c}.TCR_Strength/TCR_max;
                        now{c}=BN(now{c}, now{c}.TCR_Strength/TCR_max, now{c}.CD28_Strength/CD28_max, Cs(6)*now{c}.IL18_sum/now{c}.TCR_Stimuli, Cs(10)*now{c}.IFNg_sum/now{c}.TCR_Stimuli, Cs(11)*now{c}.TGFb_sum/now{c}.TCR_Stimuli, Cs(9)*now{c}.IL27_sum/now{c}.TCR_Stimuli, Cs(2)*now{c}.IL4_sum/now{c}.TCR_Stimuli, Cs(3)*now{c}.IL6_sum/now{c}.TCR_Stimuli, Cs(4)*now{c}.IL12_sum/now{c}.TCR_Stimuli, Cs(8)*now{c}.IL23_sum/now{c}.TCR_Stimuli, Cs(1)*now{c}.IL2_sum/now{c}.TCR_Stimuli);
                        if now{c}.Phenotype==0
                            dummy_model=Th0;
                            [MN0, Th0_MN_input, Th0_MN_output]=MN(now{c}, dummy_model, 0, protein_max_vector(1), solution0_1, solution0_2, Th0_aa_ex, Th0_aa_tran, Th0_glu_ex, Th0_glu_tran, Th0_glutaminolysis, Th0_glycolysis, Th0_lip_ex, Th0_lip_ox, Th0_lip_syn, Th0_ox_phos, Th0_pyruvate_into_mito, Th0_ATP_fluxes, Th0_ATP_coeff, Th0_AMP_fluxes, Th0_AMP_coeff, Th0_MN_input, Th0_MN_output);
                            now{c}.G1_BM_DNA=MN0(1);
                            now{c}.G1_DNA=MN0(2);
                            now{c}.G1_ATP=MN0(3);
                            now{c}.G1_AMP=MN0(4);
                            now{c}.S_BM_DNA=MN0(5);
                            now{c}.S_DNA=MN0(6);
                            now{c}.S_ATP=MN0(7);
                            now{c}.S_AMP=MN0(8);
                            now{c}.Ribosome=MN0(9);
                        elseif now{c}.Phenotype==12
                            dummy_model=Th1;
                            [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                            now{c}.G1_BM_DNA=MN1(1);
                            now{c}.G1_DNA=MN1(2);
                            now{c}.G1_ATP=MN1(3);
                            now{c}.G1_AMP=MN1(4);
                            now{c}.S_BM_DNA=MN1(5);
                            now{c}.S_DNA=MN1(6);
                            now{c}.S_ATP=MN1(7);
                            now{c}.S_AMP=MN1(8);
                            now{c}.Ribosome=MN1(9);
                        elseif now{c}.Phenotype==13
                            dummy_model=Th2;
                            [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                            now{c}.G1_BM_DNA=MN2(1);
                            now{c}.G1_DNA=MN2(2);
                            now{c}.G1_ATP=MN2(3);
                            now{c}.G1_AMP=MN2(4);
                            now{c}.S_BM_DNA=MN2(5);
                            now{c}.S_DNA=MN2(6);
                            now{c}.S_ATP=MN2(7);
                            now{c}.S_AMP=MN2(8);
                            now{c}.Ribosome=MN2(9);
                        elseif now{c}.Phenotype==14
                            dummy_model=Th17;
                            [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                            now{c}.G1_BM_DNA=MN17(1);
                            now{c}.G1_DNA=MN17(2);
                            now{c}.G1_ATP=MN17(3);
                            now{c}.G1_AMP=MN17(4);
                            now{c}.S_BM_DNA=MN17(5);
                            now{c}.S_DNA=MN17(6);
                            now{c}.S_ATP=MN17(7);
                            now{c}.S_AMP=MN17(8);
                            now{c}.Ribosome=MN17(9);
                        elseif now{c}.Phenotype==15
                            dummy_model=Treg;
                            [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);
                            now{c}.G1_BM_DNA=MNreg(1);
                            now{c}.G1_DNA=MNreg(2);
                            now{c}.G1_ATP=MNreg(3);
                            now{c}.G1_AMP=MNreg(4);
                            now{c}.S_BM_DNA=MNreg(5);
                            now{c}.S_DNA=MNreg(6);
                            now{c}.S_ATP=MNreg(7);
                            now{c}.S_AMP=MNreg(8);
                            now{c}.Ribosome=MNreg(9);
                        elseif now{c}.Phenotype==1
                            dummy_model=Th1;
                            [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                            dummy_model=Th2;
                            [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                            dummy_model=Th17;
                            [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                            dummy_model=Treg;
                            [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                            MNsum=(MN1+MN2+MN17+MNreg)/4;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);                            
                        elseif now{c}.Phenotype==2                         
                            dummy_model=Th2;
                            [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                            dummy_model=Th17;
                            [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                            dummy_model=Treg;
                            [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                            MNsum=(MN2+MN17+MNreg)/3;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);
                        elseif now{c}.Phenotype==3
                            dummy_model=Th1;
                            [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                            dummy_model=Th17;
                            [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                            dummy_model=Treg;
                            [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                            MNsum=(MN1+MN17+MNreg)/3;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);
                        elseif now{c}.Phenotype==4
                            dummy_model=Th1;
                            [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                            dummy_model=Th2;
                            [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                            dummy_model=Treg;
                            [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                            MNsum=(MN1+MN2+MNreg)/3;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);
                        elseif now{c}.Phenotype==5
                            dummy_model=Th1;
                            [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                            dummy_model=Th2;
                            [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                            dummy_model=Th17;
                            [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                            MNsum=(MN1+MN2+MN17)/3;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);    
                        elseif now{c}.Phenotype==6
                            dummy_model=Th17;
                            [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                            dummy_model=Treg;
                            [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                            MNsum=(MN17+MNreg)/2;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);
                        elseif now{c}.Phenotype==7
                            dummy_model=Th1;
                            [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                            dummy_model=Th2;
                            [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                            MNsum=(MN1+MN2)/2;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);
                        elseif now{c}.Phenotype==8
                            dummy_model=Th1;
                            [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                            dummy_model=Th17;
                            [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                            MNsum=(MN1+MN17)/2;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);    
                        elseif now{c}.Phenotype==9
                            dummy_model=Th1;
                            [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                            dummy_model=Treg;
                            [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                            MNsum=(MN1+MNreg)/2;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);    
                        elseif now{c}.Phenotype==10                     
                            dummy_model=Th2;
                            [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                            dummy_model=Treg;
                            [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                            MNsum=(MN2+MNreg)/2;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);    
                        elseif now{c}.Phenotype==11
                            dummy_model=Th2;
                            [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                            dummy_model=Th17;
                            [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                            MNsum=(MN2+MN17)/2;
                            now{c}.G1_BM_DNA=MNsum(1);
                            now{c}.G1_DNA=MNsum(2);
                            now{c}.G1_ATP=MNsum(3);
                            now{c}.G1_AMP=MNsum(4);
                            now{c}.S_BM_DNA=MNsum(5);
                            now{c}.S_DNA=MNsum(6);
                            now{c}.S_ATP=MNsum(7);
                            now{c}.S_AMP=MNsum(8);
                            now{c}.Ribosome=MNsum(9);
                        end
                        now{c}.k_tot=0;
                        now{c}.k_toln=1;
                        now{c}.k_toc=0;
                    end
                end

            %If the cell is an effector, it behaves differently depending on whether it is restimulated.
            %Concepts of Activated T Cell Death (Brenner 2008).
            %Regulation of CD4+ T-cell Contraction During Pathogen Challenge (McKinstry 2010).
            elseif now{c}.State==2
                
                %First, we need to update the migration rules. If the
                %effector cell has just arrived in the target organ, it
                %should stay there. If it has just exited the lymph node,
                %it cannot go back and will slowly move towards the target
                %organ.
                if now{c}.Location==1 %The effector cells stay in T (Alyssa's document Q3/Q5/Q6).
                    now{c}.k_tot=1;
                    now{c}.k_toln=0;
                    now{c}.k_toc=0;
                elseif now{c}.Location==3 %Update the migration probabilities to mimic its travel via the circulation system.
                    now{c}.k_tot=k_effector_CtoT;
                    now{c}.k_toln=0;
                    now{c}.k_toc=1-k_effector_CtoT;
                end
                %Second, the cell may transdifferentiate, depending on its plasticity.
                %With each division, the effector cell becomes more locked into its phenotype.
                %Effector and Memory T-Cell Differentiation: Implications for Vaccine Development (Kaech et. al. 2002).
                %Effector cells in the lungs are more activated while effector cells in the lymph node are less differentiated and produce IL2 only. Also, effector cells only exit their lymph node after several rounds of divisions.
                %CD4 T Cell Responses to Influenza Infection (Brown et. al. 2004).
                %During the first several rounds of divisions under polarising conditions, CD4+ T cells are heterogeneous, have low frequencies of cytokine production, and are phenotypically reversible.
                %Stable patterns of gene expression occur with an increasing number of cell divisions.
                %Functional and Phenotypic Plasticity of CD4+ T Cell Subsets (Caza and Landas 2015).
                %Plasticity of Human CD4+ T Cell Subsets (Geginat et. al. 2014).
                %Early Transcription and Silencing of Cytokine Genes Underlie Polarization of T Helper Cell Subsets (Grogan et. al. 2001).
                %Determine if it transdifferentiates and if it does, update its phenotype.
                %Update its metabolism by running the correct MN model.
                now{c}=BN(now{c}, now{c}.TCR_Strength/TCR_max, now{c}.CD28_Strength/CD28_max, Cs(6)*IL18(now{c}.Location), Cs(10)*IFNg(now{c}.Location), Cs(11)*TGFb(now{c}.Location), Cs(9)*IL27(now{c}.Location), Cs(2)*IL4(now{c}.Location), Cs(3)*IL6(now{c}.Location), Cs(4)*IL12(now{c}.Location), Cs(8)*IL23(now{c}.Location), Cs(1)*IL2(now{c}.Location));
                if now{c}.Phenotype==0
                    dummy_model=Th0;
                    [MN0, Th0_MN_input, Th0_MN_output]=MN(now{c}, dummy_model, 0, protein_max_vector(1), solution0_1, solution0_2, Th0_aa_ex, Th0_aa_tran, Th0_glu_ex, Th0_glu_tran, Th0_glutaminolysis, Th0_glycolysis, Th0_lip_ex, Th0_lip_ox, Th0_lip_syn, Th0_ox_phos, Th0_pyruvate_into_mito, Th0_ATP_fluxes, Th0_ATP_coeff, Th0_AMP_fluxes, Th0_AMP_coeff, Th0_MN_input, Th0_MN_output);
                    now{c}.G1_BM_DNA=MN0(1);
                    now{c}.G1_DNA=MN0(2);
                    now{c}.G1_ATP=MN0(3);
                    now{c}.G1_AMP=MN0(4);
                    now{c}.S_BM_DNA=MN0(5);
                    now{c}.S_DNA=MN0(6);
                    now{c}.S_ATP=MN0(7);
                    now{c}.S_AMP=MN0(8);
                    now{c}.Ribosome=MN0(9);
                elseif now{c}.Phenotype==12
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                    now{c}.G1_BM_DNA=MN1(1);
                    now{c}.G1_DNA=MN1(2);
                    now{c}.G1_ATP=MN1(3);
                    now{c}.G1_AMP=MN1(4);
                    now{c}.S_BM_DNA=MN1(5);
                    now{c}.S_DNA=MN1(6);
                    now{c}.S_ATP=MN1(7);
                    now{c}.S_AMP=MN1(8);
                    now{c}.Ribosome=MN1(9);
                elseif now{c}.Phenotype==13
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    now{c}.G1_BM_DNA=MN2(1);
                    now{c}.G1_DNA=MN2(2);
                    now{c}.G1_ATP=MN2(3);
                    now{c}.G1_AMP=MN2(4);
                    now{c}.S_BM_DNA=MN2(5);
                    now{c}.S_DNA=MN2(6);
                    now{c}.S_ATP=MN2(7);
                    now{c}.S_AMP=MN2(8);
                    now{c}.Ribosome=MN2(9);
                elseif now{c}.Phenotype==14
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    now{c}.G1_BM_DNA=MN17(1);
                    now{c}.G1_DNA=MN17(2);
                    now{c}.G1_ATP=MN17(3);
                    now{c}.G1_AMP=MN17(4);
                    now{c}.S_BM_DNA=MN17(5);
                    now{c}.S_DNA=MN17(6);
                    now{c}.S_ATP=MN17(7);
                    now{c}.S_AMP=MN17(8);
                    now{c}.Ribosome=MN17(9);
                elseif now{c}.Phenotype==15
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);
                    now{c}.G1_BM_DNA=MNreg(1);
                    now{c}.G1_DNA=MNreg(2);
                    now{c}.G1_ATP=MNreg(3);
                    now{c}.G1_AMP=MNreg(4);
                    now{c}.S_BM_DNA=MNreg(5);
                    now{c}.S_DNA=MNreg(6);
                    now{c}.S_ATP=MNreg(7);
                    now{c}.S_AMP=MNreg(8);
                    now{c}.Ribosome=MNreg(9);
                elseif now{c}.Phenotype==1
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN1+MN2+MN17+MNreg)/4;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);                            
                elseif now{c}.Phenotype==2                         
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN2+MN17+MNreg)/3;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==3
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN1+MN17+MNreg)/3;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==4
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN1+MN2+MNreg)/3;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==5
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    MNsum=(MN1+MN2+MN17)/3;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);    
                elseif now{c}.Phenotype==6
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN17+MNreg)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==7
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    MNsum=(MN1+MN2)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==8
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    MNsum=(MN1+MN17)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);    
                elseif now{c}.Phenotype==9
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN1+MNreg)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);    
                elseif now{c}.Phenotype==10                     
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN2+MNreg)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);    
                elseif now{c}.Phenotype==11
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    MNsum=(MN2+MN17)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                end
                
                %Third, let the cell produce biomass and switch phase if it
                %is ready. If it is already in phase 3 (G2), skip this step.
                %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                if now{c}.Phase~=3
                    now{c}=metab(now{c});
                end

                %If the cell is in G2 and it has not divided beyond its
                %limit, it may divide and return to G1.
                %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                if now{c}.Phase==3 & rand<now{c}.Pro_Cycle & now{c}.Divisions<now{c}.div_lim
                    now{c}.Divisions=now{c}.Divisions+1;
                    now{c}.Resist=1-now{c}.Divisions/now{c}.div_lim; %After a division, it loses some resistance to AICD. Concepts of Activated T Cell Death (Brenner 2008).
                    %After a division, it loses half of its material and returns to G1 (Phase 1).
                    now{c}.BM_DNA=now{c}.BM_DNA/2;
                    now{c}.DNA=now{c}.DNA/2;
                    now{c}.AMP=now{c}.AMP/2;
                    now{c}.ATP=now{c}.ATP/2;
                    now{c}.Phase=1;
                    %After a division, it loses some plasticity.
                    %However, when the division count is beyond mig_thres but below div_lim, its plasticity is zero.
                    %We know that effector cells in the lungs are more activated while effector cells in the lymph node are less differentiated and produce IL2 only. Also, effector cells only exit their lymph node after several rounds of divisions.
                    %As a result, we link the migration threshold and the transdifferentiation threshold, i.e., when the cell exits LN, it is no longer plastic.
                    %CD4 T Cell Responses to Influenza Infection (Brown et. al. 2004).
                    if now{c}.Divisions<now{c}.mig_thres
                        now{c}.Plasticity=1-now{c}.Divisions/now{c}.mig_thres;
                    else
                        now{c}.Plasticity=0;
                    end
                    %Run the BN model again to update the cell's internal state after this round of division because its resistance to AICD has decreased.
                    now{c}=BN(now{c}, now{c}.TCR_Strength/TCR_max, now{c}.CD28_Strength/CD28_max, Cs(6)*IL18(now{c}.Location), Cs(10)*IFNg(now{c}.Location), Cs(11)*TGFb(now{c}.Location), Cs(9)*IL27(now{c}.Location), Cs(2)*IL4(now{c}.Location), Cs(3)*IL6(now{c}.Location), Cs(4)*IL12(now{c}.Location), Cs(8)*IL23(now{c}.Location), Cs(1)*IL2(now{c}.Location));
                    %Run the correct MN model to update its metabolism.
                    if now{c}.Phenotype==0
                        dummy_model=Th0;
                        [MN0, Th0_MN_input, Th0_MN_output]=MN(now{c}, dummy_model, 0, protein_max_vector(1), solution0_1, solution0_2, Th0_aa_ex, Th0_aa_tran, Th0_glu_ex, Th0_glu_tran, Th0_glutaminolysis, Th0_glycolysis, Th0_lip_ex, Th0_lip_ox, Th0_lip_syn, Th0_ox_phos, Th0_pyruvate_into_mito, Th0_ATP_fluxes, Th0_ATP_coeff, Th0_AMP_fluxes, Th0_AMP_coeff, Th0_MN_input, Th0_MN_output);
                        now{c}.G1_BM_DNA=MN0(1);
                        now{c}.G1_DNA=MN0(2);
                        now{c}.G1_ATP=MN0(3);
                        now{c}.G1_AMP=MN0(4);
                        now{c}.S_BM_DNA=MN0(5);
                        now{c}.S_DNA=MN0(6);
                        now{c}.S_ATP=MN0(7);
                        now{c}.S_AMP=MN0(8);
                        now{c}.Ribosome=MN0(9);
                    elseif now{c}.Phenotype==12
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                        now{c}.G1_BM_DNA=MN1(1);
                        now{c}.G1_DNA=MN1(2);
                        now{c}.G1_ATP=MN1(3);
                        now{c}.G1_AMP=MN1(4);
                        now{c}.S_BM_DNA=MN1(5);
                        now{c}.S_DNA=MN1(6);
                        now{c}.S_ATP=MN1(7);
                        now{c}.S_AMP=MN1(8);
                        now{c}.Ribosome=MN1(9);
                    elseif now{c}.Phenotype==13
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        now{c}.G1_BM_DNA=MN2(1);
                        now{c}.G1_DNA=MN2(2);
                        now{c}.G1_ATP=MN2(3);
                        now{c}.G1_AMP=MN2(4);
                        now{c}.S_BM_DNA=MN2(5);
                        now{c}.S_DNA=MN2(6);
                        now{c}.S_ATP=MN2(7);
                        now{c}.S_AMP=MN2(8);
                        now{c}.Ribosome=MN2(9);
                    elseif now{c}.Phenotype==14
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        now{c}.G1_BM_DNA=MN17(1);
                        now{c}.G1_DNA=MN17(2);
                        now{c}.G1_ATP=MN17(3);
                        now{c}.G1_AMP=MN17(4);
                        now{c}.S_BM_DNA=MN17(5);
                        now{c}.S_DNA=MN17(6);
                        now{c}.S_ATP=MN17(7);
                        now{c}.S_AMP=MN17(8);
                        now{c}.Ribosome=MN17(9);
                    elseif now{c}.Phenotype==15
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);
                        now{c}.G1_BM_DNA=MNreg(1);
                        now{c}.G1_DNA=MNreg(2);
                        now{c}.G1_ATP=MNreg(3);
                        now{c}.G1_AMP=MNreg(4);
                        now{c}.S_BM_DNA=MNreg(5);
                        now{c}.S_DNA=MNreg(6);
                        now{c}.S_ATP=MNreg(7);
                        now{c}.S_AMP=MNreg(8);
                        now{c}.Ribosome=MNreg(9);
                    elseif now{c}.Phenotype==1
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN2+MN17+MNreg)/4;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);                            
                    elseif now{c}.Phenotype==2                         
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN2+MN17+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==3
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN17+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==4
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN2+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==5
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN1+MN2+MN17)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==6
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN17+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==7
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        MNsum=(MN1+MN2)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==8
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN1+MN17)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==9
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==10                     
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN2+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==11
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN2+MN17)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    end
                    %After division in the lymph node, if the
                    %division count has gone past mig_thres, it becomes
                    %mobile and is forced to leave for the circulation
                    %system.
                    if now{c}.Location==2 & now{c}.Divisions==now{c}.mig_thres
                        now{c}.k_tot=0;
                        now{c}.k_toln=0;
                        now{c}.k_toc=1; %If the cell has divided more than a certain number of times, it enters the circulation system. CD4 Effector T Cell Subsets in the Response to Influenza (Roman 2002).
                    end
                    %Make a copy of the cell, place it at the end of
                    %the now array, and give it the next available cell ID.
                    now(end+1)={now{c}};
                    ID=ID+1;
                    now{end}.ID=ID;
                end
                
                %With restimulation
                if rand<TCR(now{c}.Location, t)
                    %Allow the cell to decide whether it undergoes AICD.
                    %There are two pathways, hence the 'or' statement.
                    if rand<now{c}.Pro_AICD1*Ease_restim | rand<now{c}.Pro_AICD2*Ease_restim
                        now{c}.State=0;
                    end
                %Without restimulation.
                else
                    %Allow the cell to decide if it becomes a memory cell.
                    %If it does, it will have a new activation threshold.
                    %Each effector cell has the potential to form memory without additional signals apart from the withdrawl of TCR signalling; it just needs to survive long enough for that to happen.
                    %CD4(+) T Cell Effectors can Become Memory Cells with High Efficiency and Without Further Division (Hu 2001).
                    if rand<now{c}.Pro_Mem*Ease_memory & rand>now{c}.Pro_Cycle
                        now{c}.State=3;
                        now{c}.act_thres=now{c}.act_thres+act_thres_mem;
                        %If the cell has divided many times and it is
                        %highly functional, it becomes a TEM and it travels
                        %between the target organ and the circulation
                        %system. Otherwise, it becomes a TCM and stays in the LN.
                        %CD4 T-cell Memory Generation and Maintenance (Gasper 2014).
                        if now{c}.Divisions<now{c}.mig_thres
                            now{c}.k_tot=0;
                            now{c}.k_toln=1;
                            now{c}.k_toc=0;
                        else
                            now{c}.k_tot=k_effector_CtoT;
                            now{c}.k_toln=0;
                            now{c}.k_toc=1-k_effector_CtoT;
                        end
                    %Allow the cell to decide if it undergoes ACAD.
                    elseif rand<now{c}.Pro_ACAD*Ease_ACAD & rand>now{c}.Pro_Cycle
                        now{c}.State=0;                    
                    %Without restimulation, the cell is allowed to rest, so its activation level drops.
                    %CD4(+) T Cell Effectors can Become Memory Cells with High Efficiency and Without Further Division (Hu 2001).
                    else
                        now{c}.TCR_Strength=now{c}.TCR_Strength-Ease_relax*(1-TCR(now{c}.Location, t));
                        now{c}.CD28_Strength=now{c}.CD28_Strength-Ease_relax*(1-CD28(now{c}.Location, t));
                    end
                end
            
            %A memory cell is stimulated. If it has been sufficiently stimulated, it reactivates.
            %CD4 T-cell memory generation and maintenance (Gasper 2014).
            elseif now{c}.State==3
                %Restimulation.
                if rand<TCR(now{c}.Location, t)
                    now{c}.TCR_Stimuli=now{c}.TCR_Stimuli+1;
                end
                %If the cell has been restimulated enough, it reactivates.
                if now{c}.TCR_Stimuli==now{c}.act_thres
                    now{c}.State=4;
                    now{c}.TCR_Strength=TCR_max; %Reactivated memory cells respond faster and stronger, so the TCR signalling is given the maxmimum value. Qualitative Differences Between Naive and Memory T Cells (Berard et. al. 2002). CD4 T-cell Memory Generation and Maintenance (Gasper 2014).
                    now{c}.CD28_Strength=CD28_max; %Reactivated memory cells respond faster and stronger, so the CD28 signalling is given the maxmimum value. Qualitative Differences Between Naive and Memory T Cells (Berard et. al. 2002). CD4 T-cell Memory Generation and Maintenance (Gasper 2014).
                    now{c}.Ribosome=now{c}.TCR_Strength/TCR_max; %Before running the BN and MN models, we do not know the protein synthesis rate, so we will link it to the TCR strength.
                    %If the cell has not divided beyond the mig_thres limit (TCM), its migration probabilities will remain the same after reactivation, i.e. staying in the LN. However, if it has (TEM), it will return to the effector values which depend on the location. There are not any TEMs in the LN.
                    %Memory CD4+ T-cell–mediated Protection Depends on Secondary Effectors that are Distinct from and Superior to Primary Effectors (Strutt, 2012).
                    if now{c}.Divisions>=now{c}.mig_thres
                        if now{c}.Location==1
                            now{c}.k_tot=1;
                            now{c}.k_toln=0;
                            now{c}.k_toc=0;                            
                        elseif now{c}.Location==3
                            now{c}.k_tot=k_effector_CtoT;
                            now{c}.k_toln=0;
                            now{c}.k_toc=1-k_effector_CtoT;
                        end
                    end
                    %Then, run the BN model and the correct MN model with
                    %the new TCR and CD28 activity levels, as well as the
                    %current cytokine concentrations. Note that this is
                    %different from the corresponding step during the
                    %activation of a naive cell. A naive cell is perfectly
                    %plastic, so it adopts all the BN outputs. Here, if the
                    %reactivated memory cell is no longer plastic, it will
                    %only adopt the BN outputs related to metabolism, and
                    %of course the MN outputs downstream. The BN outputs
                    %related to cytokine production and transcription
                    %factors will be ignored.
                    now{c}=BN(now{c}, now{c}.TCR_Strength/TCR_max, now{c}.CD28_Strength/CD28_max, Cs(6)*IL18(now{c}.Location), Cs(10)*IFNg(now{c}.Location), Cs(11)*TGFb(now{c}.Location), Cs(9)*IL27(now{c}.Location), Cs(2)*IL4(now{c}.Location), Cs(3)*IL6(now{c}.Location), Cs(4)*IL12(now{c}.Location), Cs(8)*IL23(now{c}.Location), Cs(1)*IL2(now{c}.Location));
                    %Run the correct MN model to update its metabolism.
                    if now{c}.Phenotype==0
                        dummy_model=Th0;
                        [MN0, Th0_MN_input, Th0_MN_output]=MN(now{c}, dummy_model, 0, protein_max_vector(1), solution0_1, solution0_2, Th0_aa_ex, Th0_aa_tran, Th0_glu_ex, Th0_glu_tran, Th0_glutaminolysis, Th0_glycolysis, Th0_lip_ex, Th0_lip_ox, Th0_lip_syn, Th0_ox_phos, Th0_pyruvate_into_mito, Th0_ATP_fluxes, Th0_ATP_coeff, Th0_AMP_fluxes, Th0_AMP_coeff, Th0_MN_input, Th0_MN_output);
                        now{c}.G1_BM_DNA=MN0(1);
                        now{c}.G1_DNA=MN0(2);
                        now{c}.G1_ATP=MN0(3);
                        now{c}.G1_AMP=MN0(4);
                        now{c}.S_BM_DNA=MN0(5);
                        now{c}.S_DNA=MN0(6);
                        now{c}.S_ATP=MN0(7);
                        now{c}.S_AMP=MN0(8);
                        now{c}.Ribosome=MN0(9);
                    elseif now{c}.Phenotype==12
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                        now{c}.G1_BM_DNA=MN1(1);
                        now{c}.G1_DNA=MN1(2);
                        now{c}.G1_ATP=MN1(3);
                        now{c}.G1_AMP=MN1(4);
                        now{c}.S_BM_DNA=MN1(5);
                        now{c}.S_DNA=MN1(6);
                        now{c}.S_ATP=MN1(7);
                        now{c}.S_AMP=MN1(8);
                        now{c}.Ribosome=MN1(9);
                    elseif now{c}.Phenotype==13
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        now{c}.G1_BM_DNA=MN2(1);
                        now{c}.G1_DNA=MN2(2);
                        now{c}.G1_ATP=MN2(3);
                        now{c}.G1_AMP=MN2(4);
                        now{c}.S_BM_DNA=MN2(5);
                        now{c}.S_DNA=MN2(6);
                        now{c}.S_ATP=MN2(7);
                        now{c}.S_AMP=MN2(8);
                        now{c}.Ribosome=MN2(9);
                    elseif now{c}.Phenotype==14
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        now{c}.G1_BM_DNA=MN17(1);
                        now{c}.G1_DNA=MN17(2);
                        now{c}.G1_ATP=MN17(3);
                        now{c}.G1_AMP=MN17(4);
                        now{c}.S_BM_DNA=MN17(5);
                        now{c}.S_DNA=MN17(6);
                        now{c}.S_ATP=MN17(7);
                        now{c}.S_AMP=MN17(8);
                        now{c}.Ribosome=MN17(9);
                    elseif now{c}.Phenotype==15
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);
                        now{c}.G1_BM_DNA=MNreg(1);
                        now{c}.G1_DNA=MNreg(2);
                        now{c}.G1_ATP=MNreg(3);
                        now{c}.G1_AMP=MNreg(4);
                        now{c}.S_BM_DNA=MNreg(5);
                        now{c}.S_DNA=MNreg(6);
                        now{c}.S_ATP=MNreg(7);
                        now{c}.S_AMP=MNreg(8);
                        now{c}.Ribosome=MNreg(9);
                    elseif now{c}.Phenotype==1
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN2+MN17+MNreg)/4;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);                            
                    elseif now{c}.Phenotype==2                         
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN2+MN17+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==3
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN17+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==4
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN2+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==5
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN1+MN2+MN17)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==6
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN17+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==7
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        MNsum=(MN1+MN2)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==8
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN1+MN17)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==9
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==10                     
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN2+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==11
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN2+MN17)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    end
                end
 
            %A reactivated memory cell acts like one derived from a naive cell.
            %Migration, Maintenance and Recall of Memory T Cells in Peripheral Tissues (Woodland, 2009).
            else
                
                %First, we need to update the migration rules. If the
                %cell has just arrived in the target organ, it should stay there.
                %If it has just exited the lymph node, it cannot go back and will slowly move towards the target organ.
                if now{c}.Location==1 %The effector cells stay in T (Alyssa's document Q3/Q5/Q6).
                    now{c}.k_tot=1;
                    now{c}.k_toln=0;
                    now{c}.k_toc=0;
                elseif now{c}.Location==3 %Update the migration probabilities to mimic its travel via the circulation system.
                    now{c}.k_tot=k_effector_CtoT;
                    now{c}.k_toln=0;
                    now{c}.k_toc=1-k_effector_CtoT;
                end
                
                %Second, the cell may transdifferentiate, depending on its plasticity.
                %TEMs are limited to the periphery, lineage-specific, produce a wide range of cytokines, and do not proliferate as well as TCMs. They are derived from cytokine-producing cells.
                %TCMs are limited to secondary lymphoid organs, produce IL2 only, and proliferate well.
                %Qualitative Differences Between Naive and Memory T Cells (Berard et. al. 2002). 
                %CD4 T-cell Memory Generation and Maintenance (Gasper 2014).
                %As a result, we equate TCMs with plastic cells, TEMs with non-plastic cells.
                %Determine if it transdifferentiates and if it does, update its phenotype.
                %Update its metabolism by running the correct MN model.
                now{c}=BN(now{c}, now{c}.TCR_Strength/TCR_max, now{c}.CD28_Strength/CD28_max, Cs(6)*IL18(now{c}.Location), Cs(10)*IFNg(now{c}.Location), Cs(11)*TGFb(now{c}.Location), Cs(9)*IL27(now{c}.Location), Cs(2)*IL4(now{c}.Location), Cs(3)*IL6(now{c}.Location), Cs(4)*IL12(now{c}.Location), Cs(8)*IL23(now{c}.Location), Cs(1)*IL2(now{c}.Location));
                if now{c}.Phenotype==0
                    dummy_model=Th0;
                    [MN0, Th0_MN_input, Th0_MN_output]=MN(now{c}, dummy_model, 0, protein_max_vector(1), solution0_1, solution0_2, Th0_aa_ex, Th0_aa_tran, Th0_glu_ex, Th0_glu_tran, Th0_glutaminolysis, Th0_glycolysis, Th0_lip_ex, Th0_lip_ox, Th0_lip_syn, Th0_ox_phos, Th0_pyruvate_into_mito, Th0_ATP_fluxes, Th0_ATP_coeff, Th0_AMP_fluxes, Th0_AMP_coeff, Th0_MN_input, Th0_MN_output);
                    now{c}.G1_BM_DNA=MN0(1);
                    now{c}.G1_DNA=MN0(2);
                    now{c}.G1_ATP=MN0(3);
                    now{c}.G1_AMP=MN0(4);
                    now{c}.S_BM_DNA=MN0(5);
                    now{c}.S_DNA=MN0(6);
                    now{c}.S_ATP=MN0(7);
                    now{c}.S_AMP=MN0(8);
                    now{c}.Ribosome=MN0(9);
                elseif now{c}.Phenotype==12
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                    now{c}.G1_BM_DNA=MN1(1);
                    now{c}.G1_DNA=MN1(2);
                    now{c}.G1_ATP=MN1(3);
                    now{c}.G1_AMP=MN1(4);
                    now{c}.S_BM_DNA=MN1(5);
                    now{c}.S_DNA=MN1(6);
                    now{c}.S_ATP=MN1(7);
                    now{c}.S_AMP=MN1(8);
                    now{c}.Ribosome=MN1(9);
                elseif now{c}.Phenotype==13
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    now{c}.G1_BM_DNA=MN2(1);
                    now{c}.G1_DNA=MN2(2);
                    now{c}.G1_ATP=MN2(3);
                    now{c}.G1_AMP=MN2(4);
                    now{c}.S_BM_DNA=MN2(5);
                    now{c}.S_DNA=MN2(6);
                    now{c}.S_ATP=MN2(7);
                    now{c}.S_AMP=MN2(8);
                    now{c}.Ribosome=MN2(9);
                elseif now{c}.Phenotype==14
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    now{c}.G1_BM_DNA=MN17(1);
                    now{c}.G1_DNA=MN17(2);
                    now{c}.G1_ATP=MN17(3);
                    now{c}.G1_AMP=MN17(4);
                    now{c}.S_BM_DNA=MN17(5);
                    now{c}.S_DNA=MN17(6);
                    now{c}.S_ATP=MN17(7);
                    now{c}.S_AMP=MN17(8);
                    now{c}.Ribosome=MN17(9);
                elseif now{c}.Phenotype==15
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);
                    now{c}.G1_BM_DNA=MNreg(1);
                    now{c}.G1_DNA=MNreg(2);
                    now{c}.G1_ATP=MNreg(3);
                    now{c}.G1_AMP=MNreg(4);
                    now{c}.S_BM_DNA=MNreg(5);
                    now{c}.S_DNA=MNreg(6);
                    now{c}.S_ATP=MNreg(7);
                    now{c}.S_AMP=MNreg(8);
                    now{c}.Ribosome=MNreg(9);
                elseif now{c}.Phenotype==1
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN1+MN2+MN17+MNreg)/4;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);                            
                elseif now{c}.Phenotype==2                         
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN2+MN17+MNreg)/3;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==3
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN1+MN17+MNreg)/3;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==4
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN1+MN2+MNreg)/3;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==5
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    MNsum=(MN1+MN2+MN17)/3;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);    
                elseif now{c}.Phenotype==6
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN17+MNreg)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==7
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    MNsum=(MN1+MN2)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                elseif now{c}.Phenotype==8
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    MNsum=(MN1+MN17)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);    
                elseif now{c}.Phenotype==9
                    dummy_model=Th1;
                    [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN1+MNreg)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);    
                elseif now{c}.Phenotype==10                     
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Treg;
                    [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                    MNsum=(MN2+MNreg)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);    
                elseif now{c}.Phenotype==11
                    dummy_model=Th2;
                    [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                    dummy_model=Th17;
                    [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                    MNsum=(MN2+MN17)/2;
                    now{c}.G1_BM_DNA=MNsum(1);
                    now{c}.G1_DNA=MNsum(2);
                    now{c}.G1_ATP=MNsum(3);
                    now{c}.G1_AMP=MNsum(4);
                    now{c}.S_BM_DNA=MNsum(5);
                    now{c}.S_DNA=MNsum(6);
                    now{c}.S_ATP=MNsum(7);
                    now{c}.S_AMP=MNsum(8);
                    now{c}.Ribosome=MNsum(9);
                end
                
                %Third, let the cell produce biomass and switch phase if it
                %is ready. If it is already in phase 3 (G2), skip this step.
                %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                if now{c}.Phase~=3
                    now{c}=metab(now{c});
                end

                %If the cell is in G2 and it has not divided beyond its
                %limit, it may divide and return to G1.
                %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                if now{c}.Phase==3 & rand<now{c}.Pro_Cycle & now{c}.Divisions<now{c}.div_lim
                    now{c}.Divisions=now{c}.Divisions+1;
                    now{c}.Resist=1-now{c}.Divisions/now{c}.div_lim; %After a division, it loses some resistance to AICD. Concepts of Activated T Cell Death (Brenner 2008).
                    %After a division, it loses half of its material and returns to G1 (Phase 1).
                    now{c}.BM_DNA=now{c}.BM_DNA/2;
                    now{c}.DNA=now{c}.DNA/2;
                    now{c}.AMP=now{c}.AMP/2;
                    now{c}.ATP=now{c}.ATP/2;
                    now{c}.Phase=1;
                    %After a division, it loses some plasticity.
                    %However, when the division count is beyond mig_thres but below div_lim, its plasticity is zero.
                    %TEMs are limited to the periphery, lineage-specific, produce a wide range of cytokines, and do not proliferate as well as TCMs. They are derived from cytokine-producing cells.
                    %TCMs are limited to secondary lymphoid organs, produce IL2 only, and proliferate well.
                    %Qualitative Differences Between Naive and Memory T Cells (Berard et. al. 2002). 
                    %CD4 T-cell Memory Generation and Maintenance (Gasper 2014).
                    %As a result, we link the migration threshold and the transdifferentiation threshold, i.e., when the cell exits LN, it is no longer plastic.                   
                    if now{c}.Divisions<now{c}.mig_thres
                        now{c}.Plasticity=1-now{c}.Divisions/now{c}.mig_thres;
                    else
                        now{c}.Plasticity=0;
                    end
                    %Run the BN model again to update the cell's internal state after this round of division because its resistance to AICD has decreased.
                    now{c}=BN(now{c}, now{c}.TCR_Strength/TCR_max, now{c}.CD28_Strength/CD28_max, Cs(6)*IL18(now{c}.Location), Cs(10)*IFNg(now{c}.Location), Cs(11)*TGFb(now{c}.Location), Cs(9)*IL27(now{c}.Location), Cs(2)*IL4(now{c}.Location), Cs(3)*IL6(now{c}.Location), Cs(4)*IL12(now{c}.Location), Cs(8)*IL23(now{c}.Location), Cs(1)*IL2(now{c}.Location));
                    %Run the correct MN model to update its metabolism.
                    if now{c}.Phenotype==0
                        dummy_model=Th0;
                        [MN0, Th0_MN_input, Th0_MN_output]=MN(now{c}, dummy_model, 0, protein_max_vector(1), solution0_1, solution0_2, Th0_aa_ex, Th0_aa_tran, Th0_glu_ex, Th0_glu_tran, Th0_glutaminolysis, Th0_glycolysis, Th0_lip_ex, Th0_lip_ox, Th0_lip_syn, Th0_ox_phos, Th0_pyruvate_into_mito, Th0_ATP_fluxes, Th0_ATP_coeff, Th0_AMP_fluxes, Th0_AMP_coeff, Th0_MN_input, Th0_MN_output);
                        now{c}.G1_BM_DNA=MN0(1);
                        now{c}.G1_DNA=MN0(2);
                        now{c}.G1_ATP=MN0(3);
                        now{c}.G1_AMP=MN0(4);
                        now{c}.S_BM_DNA=MN0(5);
                        now{c}.S_DNA=MN0(6);
                        now{c}.S_ATP=MN0(7);
                        now{c}.S_AMP=MN0(8);
                        now{c}.Ribosome=MN0(9);
                    elseif now{c}.Phenotype==12
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                        now{c}.G1_BM_DNA=MN1(1);
                        now{c}.G1_DNA=MN1(2);
                        now{c}.G1_ATP=MN1(3);
                        now{c}.G1_AMP=MN1(4);
                        now{c}.S_BM_DNA=MN1(5);
                        now{c}.S_DNA=MN1(6);
                        now{c}.S_ATP=MN1(7);
                        now{c}.S_AMP=MN1(8);
                        now{c}.Ribosome=MN1(9);
                    elseif now{c}.Phenotype==13
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        now{c}.G1_BM_DNA=MN2(1);
                        now{c}.G1_DNA=MN2(2);
                        now{c}.G1_ATP=MN2(3);
                        now{c}.G1_AMP=MN2(4);
                        now{c}.S_BM_DNA=MN2(5);
                        now{c}.S_DNA=MN2(6);
                        now{c}.S_ATP=MN2(7);
                        now{c}.S_AMP=MN2(8);
                        now{c}.Ribosome=MN2(9);
                    elseif now{c}.Phenotype==14
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        now{c}.G1_BM_DNA=MN17(1);
                        now{c}.G1_DNA=MN17(2);
                        now{c}.G1_ATP=MN17(3);
                        now{c}.G1_AMP=MN17(4);
                        now{c}.S_BM_DNA=MN17(5);
                        now{c}.S_DNA=MN17(6);
                        now{c}.S_ATP=MN17(7);
                        now{c}.S_AMP=MN17(8);
                        now{c}.Ribosome=MN17(9);
                    elseif now{c}.Phenotype==15
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);
                        now{c}.G1_BM_DNA=MNreg(1);
                        now{c}.G1_DNA=MNreg(2);
                        now{c}.G1_ATP=MNreg(3);
                        now{c}.G1_AMP=MNreg(4);
                        now{c}.S_BM_DNA=MNreg(5);
                        now{c}.S_DNA=MNreg(6);
                        now{c}.S_ATP=MNreg(7);
                        now{c}.S_AMP=MNreg(8);
                        now{c}.Ribosome=MNreg(9);
                    elseif now{c}.Phenotype==1
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN2+MN17+MNreg)/4;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);                            
                    elseif now{c}.Phenotype==2                         
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN2+MN17+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==3
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN17+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==4
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MN2+MNreg)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==5
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN1+MN2+MN17)/3;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==6
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN17+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==7
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        MNsum=(MN1+MN2)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    elseif now{c}.Phenotype==8
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN1+MN17)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==9
                        dummy_model=Th1;
                        [MN1, Th1_MN_input, Th1_MN_output]=MN(now{c}, dummy_model, 12, protein_max_vector(2), solution12_1, solution12_2, Th1_aa_ex, Th1_aa_tran, Th1_glu_ex, Th1_glu_tran, Th1_glutaminolysis, Th1_glycolysis, Th1_lip_ex, Th1_lip_ox, Th1_lip_syn, Th1_ox_phos, Th1_pyruvate_into_mito, Th1_ATP_fluxes, Th1_ATP_coeff, Th1_AMP_fluxes, Th1_AMP_coeff, Th1_MN_input, Th1_MN_output);                            
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN1+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==10                     
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Treg;
                        [MNreg, Treg_MN_input, Treg_MN_output]=MN(now{c}, dummy_model, 15, protein_max_vector(5), solution15_1, solution15_2, Treg_aa_ex, Treg_aa_tran, Treg_glu_ex, Treg_glu_tran, Treg_glutaminolysis, Treg_glycolysis, Treg_lip_ex, Treg_lip_ox, Treg_lip_syn, Treg_ox_phos, Treg_pyruvate_into_mito, Treg_ATP_fluxes, Treg_ATP_coeff, Treg_AMP_fluxes, Treg_AMP_coeff, Treg_MN_input, Treg_MN_output);                            
                        MNsum=(MN2+MNreg)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);    
                    elseif now{c}.Phenotype==11
                        dummy_model=Th2;
                        [MN2, Th2_MN_input, Th2_MN_output]=MN(now{c}, dummy_model, 13, protein_max_vector(3), solution13_1, solution13_2, Th2_aa_ex, Th2_aa_tran, Th2_glu_ex, Th2_glu_tran, Th2_glutaminolysis, Th2_glycolysis, Th2_lip_ex, Th2_lip_ox, Th2_lip_syn, Th2_ox_phos, Th2_pyruvate_into_mito, Th2_ATP_fluxes, Th2_ATP_coeff, Th2_AMP_fluxes, Th2_AMP_coeff, Th2_MN_input, Th2_MN_output);
                        dummy_model=Th17;
                        [MN17, Th17_MN_input, Th17_MN_output]=MN(now{c}, dummy_model, 14, protein_max_vector(4), solution14_1, solution14_2, Th17_aa_ex, Th17_aa_tran, Th17_glu_ex, Th17_glu_tran, Th17_glutaminolysis, Th17_glycolysis, Th17_lip_ex, Th17_lip_ox, Th17_lip_syn, Th17_ox_phos, Th17_pyruvate_into_mito, Th17_ATP_fluxes, Th17_ATP_coeff, Th17_AMP_fluxes, Th17_AMP_coeff, Th17_MN_input, Th17_MN_output);
                        MNsum=(MN2+MN17)/2;
                        now{c}.G1_BM_DNA=MNsum(1);
                        now{c}.G1_DNA=MNsum(2);
                        now{c}.G1_ATP=MNsum(3);
                        now{c}.G1_AMP=MNsum(4);
                        now{c}.S_BM_DNA=MNsum(5);
                        now{c}.S_DNA=MNsum(6);
                        now{c}.S_ATP=MNsum(7);
                        now{c}.S_AMP=MNsum(8);
                        now{c}.Ribosome=MNsum(9);
                    end
                    %After division in the lymph node, if the
                    %division count has gone past mig_thres, it becomes
                    %mobile and is forced to leave for the circulation
                    %system.
                    if now{c}.Location==2 & now{c}.Divisions==now{c}.mig_thres
                        now{c}.k_tot=0;
                        now{c}.k_toln=0;
                        now{c}.k_toc=1; %If the cell has divided more than a certain number of times, it enters the circulation system. CD4 Effector T Cell Subsets in the Response to Influenza (Roman 2002).
                    end
                    %Make a copy of the cell, place it at the end of
                    %the now array, and give it the next available cell ID.
                    now(end+1)={now{c}};
                    ID=ID+1;
                    now{end}.ID=ID;
                end
                
                %With restimulation
                if rand<TCR(now{c}.Location, t)
                    %Allow the cell to decide whether it undergoes AICD.
                    %There are two pathways, hence the 'or' statement.
                    if rand<now{c}.Pro_AICD1*Ease_restim | rand<now{c}.Pro_AICD2*Ease_restim
                        now{c}.State=0;
                    end
                %Without restimulation.
                else
                    %Allow the cell to decide if it becomes a memory cell.
                    %If it does, it will have a new activation threshold.
                    %Each effector cell has the potential to form memory without additional signals apart from the withdrawl of TCR signalling; it just needs to survive long enough for that to happen.
                    %CD4(+) T Cell Effectors can Become Memory Cells with High Efficiency and Without Further Division (Hu 2001).
                    %Secondary challenges shift the balance from TCMs to TEMs.
                    %CD4 T-cell Memory Generation and Maintenance (Gasper 2014). 
                    if rand<now{c}.Pro_Mem*Ease_memory & rand>now{c}.Pro_Cycle
                        now{c}.State=3;
                        now{c}.act_thres=now{c}.act_thres+act_thres_mem;
                        %If the cell has divided many times and it is
                        %highly functional, it becomes a TEM and it travels
                        %between the target organ and the circulation
                        %system. Otherwise, it becomes a TCM and stays in the LN.
                        %CD4 T-cell Memory Generation and Maintenance (Gasper 2014).                        
                        if now{c}.Divisions<now{c}.mig_thres
                            now{c}.k_tot=0;
                            now{c}.k_toln=1;
                            now{c}.k_toc=0;
                        else
                            now{c}.k_tot=k_effector_CtoT;
                            now{c}.k_toln=0;
                            now{c}.k_toc=1-k_effector_CtoT;
                        end
                    %Allow the cell to decide if it undergoes ACAD.
                    elseif rand<now{c}.Pro_ACAD*Ease_ACAD & rand>now{c}.Pro_Cycle
                        now{c}.State=0;                    
                    %Without restimulation, the cell is allowed to rest, so its activation level drops.
                    %CD4(+) T Cell Effectors can Become Memory Cells with High Efficiency and Without Further Division (Hu 2001).
                    else
                        now{c}.TCR_Strength=now{c}.TCR_Strength-Ease_relax*(1-TCR(now{c}.Location, t));
                        now{c}.CD28_Strength=now{c}.CD28_Strength-Ease_relax*(1-CD28(now{c}.Location, t));
                    end
                end

            end
            
        end
        
        %Count the number of living cells at the end of the time step. Also, record their positions in the now array. Redistribute each living cell in space. 
        %If there are living cells in the now array at the end of the time step, just remove the dead ones from the now array.
        %If there are only dead cells in the now array, delete the whole array, create a new one, and populate it with the initial population of naive cells.
        %Note that the cells created during this time step are also considered.
        living_count=0;
        living_position=[];
        for c=1:size(now, 1)
            if now{c}.State~=0
                living_count=living_count+1;
                living_position=[living_position, c];
                dummy=rand;
                if dummy<now{c}.k_tot
                    now{c}.Location=1;
                elseif dummy<now{c}.k_tot+now{c}.k_toln
                    now{c}.Location=2;
                else
                    now{c}.Location=3;
                end                 
            end
        end
        if living_count~=0
            now=now(living_position);
        else
            clear now;
            now=cell(Nt0+Nln0+Nc0, 1);
            for c=1:Nt0
                ID=ID+1;
                now{c}=CD4(ID, 1, step_size, 1-Nln0/Nt0, Nln0/Nt0/2, Nln0/Nt0/2);
            end
            for c=Nt0+1:Nt0+Nln0
                ID=ID+1;
                now{c}=CD4(ID, 2, step_size, 0, 0.5, 0.5);
            end
            for c=Nt0+Nln0+1:Nt0+Nln0+Nc0
                ID=ID+1;
                now{c}=CD4(ID, 3, step_size, Nln0/Nc0, 0, 1-Nln0/Nc0);
            end
        end
        
        %Count the number of naive cells in each compartment. If the counts
        %in all three compartments are below the corresponding initial
        %counts, refill the compartments.
        naive_count=zeros(1, 3);
        for c=1:size(now, 1)
            if now{c}.State==1
                naive_count(1, now{c}.Location)=naive_count(1, now{c}.Location)+1;
            end
        end
        Nt0_diff=Nt0-naive_count(1, 1);
        Nln0_diff=Nln0-naive_count(1, 2);
        Nc0_diff=Nc0-naive_count(1, 3);
        if Nt0_diff>0 & Nln0_diff>0 & Nc0_diff>0 & (naive_count(1, 1)+naive_count(1, 2)+naive_count(1, 3))<(Nt0+Nln0+Nc0)*Refill_thres
            now_refill=cell(Nt0_diff+Nln0_diff+Nc0_diff, 1);
            for c=1:Nt0_diff
                ID=ID+1;
                now_refill{c}=CD4(ID, 1, step_size, 1-Nln0/Nt0, Nln0/Nt0/2, Nln0/Nt0/2);
            end
            for c=Nt0_diff+1:Nt0_diff+Nln0_diff
                ID=ID+1;
                now_refill{c}=CD4(ID, 2, step_size, 0, 0.5, 0.5);
            end
            for c=Nt0_diff+Nln0_diff+1:Nt0_diff+Nln0_diff+Nc0_diff
                ID=ID+1;
                now_refill{c}=CD4(ID, 3, step_size, Nln0/Nc0, 0, 1-Nln0/Nc0);
            end
            now=[now;now_refill];
        end

        %Calculate the cytokine concentrations at the end of the time step. The
        %inputs are all from the beginning of the time step. This task must
        %take place before we calculate how many cells can produce each cytokine 
        %in each compartment at the end of the time step.
        [IL2(1), IL2(2), IL2(3)]=cytokines(IL20(1), IL20(2), IL20(3), IL2(1), IL2(2), IL2(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), N_IL2(1), N_IL2(2), N_IL2(3), tau, Cs(1), kdeg, Qa, Qb, Vt, Vln, Vc, PIL2(1), PIL2in(1), PIL2ex(1), PIL2(2), PIL2in(2), PIL2ex(2), PIL2(3), PIL2in(3), PIL2ex(3));
        [IL4(1), IL4(2), IL4(3)]=cytokines(IL40(1), IL40(2), IL40(3), IL4(1), IL4(2), IL4(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), N_IL4(1), N_IL4(2), N_IL4(3), tau, Cs(2), kdeg, Qa, Qb, Vt, Vln, Vc, PIL4(1), PIL4in(1), PIL4ex(1), PIL4(2), PIL4in(2), PIL4ex(2), PIL4(3), PIL4in(3), PIL4ex(3));
        [IL6(1), IL6(2), IL6(3)]=cytokines(IL60(1), IL60(2), IL60(3), IL6(1), IL6(2), IL6(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), N_IL6(1), N_IL6(2), N_IL6(3), tau, Cs(3), kdeg, Qa, Qb, Vt, Vln, Vc, PIL6(1), PIL6in(1), PIL6ex(1), PIL6(2), PIL6in(2), PIL6ex(2), PIL6(3), PIL6in(3), PIL6ex(3));
        [IL12(1), IL12(2), IL12(3)]=cytokines(IL120(1), IL120(2), IL120(3), IL12(1), IL12(2), IL12(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), 0, 0, 0, tau, Cs(4), kdeg, Qa, Qb, Vt, Vln, Vc, PIL12(1), PIL12in(1), PIL12ex(1), PIL12(2), PIL12in(2), PIL12ex(2), PIL12(3), PIL12in(3), PIL12ex(3));
        [IL17(1), IL17(2), IL17(3)]=cytokines(IL170(1), IL170(2), IL170(3), IL17(1), IL17(2), IL17(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), N_IL17(1), N_IL17(2), N_IL17(3), tau, Cs(5), kdeg, Qa, Qb, Vt, Vln, Vc, PIL17(1), PIL17in(1), PIL17ex(1), PIL17(2), PIL17in(2), PIL17ex(2), PIL17(3), PIL17in(3), PIL17ex(3));
        [IL18(1), IL18(2), IL18(3)]=cytokines(IL180(1), IL180(2), IL180(3), IL18(1), IL18(2), IL18(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), 0, 0, 0, tau, Cs(6), kdeg, Qa, Qb, Vt, Vln, Vc, PIL18(1), PIL18in(1), PIL18ex(1), PIL18(2), PIL18in(2), PIL18ex(2), PIL18(3), PIL18in(3), PIL18ex(3));
        [IL21(1), IL21(2), IL21(3)]=cytokines(IL210(1), IL210(2), IL210(3), IL21(1), IL21(2), IL21(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), N_IL21(1), N_IL21(2), N_IL21(3), tau, Cs(7), kdeg, Qa, Qb, Vt, Vln, Vc, PIL21(1), PIL21in(1), PIL21ex(1), PIL21(2), PIL21in(2), PIL21ex(2), PIL21(3), PIL21in(3), PIL21ex(3));
        [IL23(1), IL23(2), IL23(3)]=cytokines(IL230(1), IL230(2), IL230(3), IL23(1), IL23(2), IL23(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), 0, 0, 0, tau, Cs(8), kdeg, Qa, Qb, Vt, Vln, Vc, PIL23(1), PIL23in(1), PIL23ex(1), PIL23(2), PIL23in(2), PIL23ex(2), PIL23(3), PIL23in(3), PIL23ex(3));
        [IL27(1), IL27(2), IL27(3)]=cytokines(IL270(1), IL270(2), IL270(3), IL27(1), IL27(2), IL27(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), 0, 0, 0, tau, Cs(9), kdeg, Qa, Qb, Vt, Vln, Vc, PIL27(1), PIL27in(1), PIL27ex(1), PIL27(2), PIL27in(2), PIL27ex(2), PIL27(3), PIL27in(3), PIL27ex(3));
        [IFNg(1), IFNg(2), IFNg(3)]=cytokines(IFNg0(1), IFNg0(2), IFNg0(3), IFNg(1), IFNg(2), IFNg(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), N_IFNg(1), N_IFNg(2), N_IFNg(3), tau, Cs(10), kdeg, Qa, Qb, Vt, Vln, Vc, PIFNg(1), PIFNgin(1), PIFNgex(1), PIFNg(2), PIFNgin(2), PIFNgex(2), PIFNg(3), PIFNgin(3), PIFNgex(3));
        [TGFb(1), TGFb(2), TGFb(3)]=cytokines(TGFb0(1), TGFb0(2), TGFb0(3), TGFb(1), TGFb(2), TGFb(3), Antigen_target(t), Antigen_LN(t), Antigen_circulation(t), 0, 0, 0, tau, Cs(11), kdeg, Qa, Qb, Vt, Vln, Vc, PTGFb(1), PTGFbin(1), PTGFbex(1), PTGFb(2), PTGFbin(2), PTGFbex(2), PTGFb(3), PTGFbin(3), PTGFbex(3));

        %Work out how many cells can produce each cytokine in each compartment at the end of the time step. These counts will be used as inputs in the next time step.
        %Only effector cells and reactivated memory cells can produce cytokines. In both cases, cytokine production increases with decreasing plasticity. 
        %During the first several rounds of divisions under polarising conditions, CD4+ T cells are heterogeneous, have low frequencies of cytokine production, and are phenotypically reversible.
        %Stable patterns of gene expression occur with an increasing number of cell divisions.
        %Functional and Phenotypic Plasticity of CD4+ T Cell Subsets (Caza and Landas 2015).
        %Plasticity of Human CD4+ T Cell Subsets (Geginat et. al. 2014).
        %Early Transcription and Silencing of Cytokine Genes Underlie Polarization of T Helper Cell Subsets (Grogan et. al. 2001).
        %TEMs are limited to the periphery, lineage-specific, produce a wide range of cytokines, and do not proliferate as well as TCMs. They are derived from cytokine-producing cells.
        %TCMs are limited to secondary lymphoid organs, produce IL2 only, and proliferate well.
        %Qualitative Differences Between Naive and Memory T Cells (Berard et. al. 2002). 
        %CD4 T-cell Memory Generation and Maintenance (Gasper 2014).
        N_IFNg=[0, 0, 0];
        N_IL4=[0, 0, 0];
        N_IL6=[0, 0, 0];
        N_IL2=[0, 0, 0];
        N_IL17=[0, 0, 0];
        N_IL21=[0, 0, 0];
        %Consider each cell.
        for c=1:size(now, 1)
            if now{c}.State==2 | now{c}.State==4
                dummy=now{c}.Location;
                dummy2=1-now{c}.Plasticity;
                %Remember to convert the activity levels of the cytokine
                %production nodes into their empirical cumulative
                %distribution function values. However, the activity levels
                %recorded inside the cell will not change.
                scaled_IL2=Cytokines_conversion(now{c}.Pro_IL2, 1);
                scaled_IL4=Cytokines_conversion(now{c}.Pro_IL4, 2);
                scaled_IL6=Cytokines_conversion(now{c}.Pro_IL6, 3);
                scaled_IL17=Cytokines_conversion(now{c}.Pro_IL17, 4);
                scaled_IL21=Cytokines_conversion(now{c}.Pro_IL21, 5);
                scaled_IFNg=Cytokines_conversion(now{c}.Pro_IFNg, 6);
                %Depending on the cell's plasticity, it either produces
                %lineage-specific cytokines or IL2 only.
                if rand<dummy2
                    N_IFNg(dummy)=N_IFNg(dummy)+(rand<scaled_IFNg);
                    N_IL4(dummy)=N_IL4(dummy)+(rand<scaled_IL4);
                    N_IL6(dummy)=N_IL6(dummy)+(rand<scaled_IL6);
                    N_IL2(dummy)=N_IL2(dummy)+(rand<scaled_IL2);
                    N_IL17(dummy)=N_IL17(dummy)+(rand<scaled_IL17);
                    N_IL21(dummy)=N_IL21(dummy)+(rand<scaled_IL21);
                else
                    N_IL2(dummy)=N_IL2(dummy)+now{c}.Plasticity;
                end
            end
        end

        %Record the end of the time step: the entire now array and the cytokine concentrations.
        %In log_cytokines, entries 1-11 are the concentrations in T; 12-22, LN; 23-33, C.
        clear log_cytokines;
        cytokines_t=[IL2(1)*Cs(1), IL4(1)*Cs(2), IL6(1)*Cs(3), IL12(1)*Cs(4), IL17(1)*Cs(5), IL18(1)*Cs(6), IL21(1)*Cs(7), IL23(1)*Cs(8), IL27(1)*Cs(9), IFNg(1)*Cs(10), TGFb(1)*Cs(11)];
        cytokines_ln=[IL2(2)*Cs(1), IL4(2)*Cs(2), IL6(2)*Cs(3), IL12(2)*Cs(4), IL17(2)*Cs(5), IL18(2)*Cs(6), IL21(2)*Cs(7), IL23(2)*Cs(8), IL27(2)*Cs(9), IFNg(2)*Cs(10), TGFb(2)*Cs(11)];
        cytokines_c=[IL2(3)*Cs(1), IL4(3)*Cs(2), IL6(3)*Cs(3), IL12(3)*Cs(4), IL17(3)*Cs(5), IL18(3)*Cs(6), IL21(3)*Cs(7), IL23(3)*Cs(8), IL27(3)*Cs(9), IFNg(3)*Cs(10), TGFb(3)*Cs(11)];
        log_cytokines=[cytokines_t, cytokines_ln, cytokines_c];
        fullname=['Simulation', num2str(name), '_', num2str(t), '.mat'];
        save(fullname, 'now', 'log_cytokines');
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the updated libraries of MN solutions.
save('Th0_MNlibrary', 'Th0_MN_input', 'Th0_MN_output');
save('Th1_MNlibrary', 'Th1_MN_input', 'Th1_MN_output');
save('Th2_MNlibrary', 'Th2_MN_input', 'Th2_MN_output');
save('Th17_MNlibrary', 'Th17_MN_input', 'Th17_MN_output');
save('Treg_MNlibrary', 'Treg_MN_input', 'Treg_MN_output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Timekeeping ends.
toc
