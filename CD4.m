classdef CD4
    %This class represents a generic CD4+ T cell.
    
    properties
        
        ID %Unique ID number of the cell.
        step_size %The size of a time step (hours).
        div_lim %Division limit.
        mig_thres %Number of divisions the cell must undergo before it can exit LN.
        act_thres %Number of stimuli the cell must accumulate before its activation or reactivation.
        State %Cell state: 1=naive, 2=effector, 3=memory, 4=reactivated memory, and 0=dead.
        BN_count %Number of times the BN model has been run for this cell.
        Phenotype %Cellular phenotype: 0=none (Th0), 1=Tbet-GATA3-RORgt-Foxp3, 2=GATA3-RORgt-Foxp3, 3=Tbet-RORgt-Foxp3, 4=Tbet-GATA3-Foxp3, 5=Tbet-GATA3-RORgt, 6=RORgt-Foxp3, 7=Tbet-GATA3, 8=Tbet-RORgt, 9=Tbet-Foxp3, 10=GATA3-Foxp3, 11=GATA3-RORgt, 12=T-bet (Th1), 13=GATA3 (Th2), 14=RORgt (Th17), and 15=Foxp3 (Treg).
        Phase %Cell cycle phase (1=G1, 2=S, and 3=G2).
        Location %Compartment: 1=target, 2=lymph node, and 3=circulation.
        k_tot %Probability of migration to T.
        k_toln %Probability of migration to LN.
        k_toc %Probability of migration to C.
        Divisions %Number of divisions the cell has undergone.
        Plasticity %Ability to transdifferentiate.
        TCR_Stimuli %Number of TCR stimuli accumulated.
        TCR_Strength %Sum of the magnitudes of the TCR stimuli accumulated.
        CD28_Strength %Sum of the magnitudes of the CD28 stimuli accumulated.
        mTORC1_t %Prolonged activity of mTORC1.
        Ribosome %Activity level of ribosome.
        Resist %Resistance to AICD.
        IL2_sum %External cytokine signals accumulated.
        IL4_sum %External cytokine signals accumulated.
        IL6_sum %External cytokine signals accumulated.
        IL12_sum %External cytokine signals accumulated.
        IL18_sum %External cytokine signals accumulated.
        IL23_sum %External cytokine signals accumulated.
        IL27_sum %External cytokine signals accumulated.
        IFNg_sum %External cytokine signals accumulated.
        TGFb_sum %External cytokine signals accumulated.
        Pro_IL6R %Probability that this internal node is on.
        Pro_Jak1 %Probability that this internal node is on.
        Pro_Gluta %Probability that this internal node is on.
        Pro_IL23R %Probability that this internal node is on.
        Pro_Akt %Probability that this internal node is on.
        Pro_cMyc %Probability that this internal node is on. 
        Pro_FasLR %Probability that this internal node is on.
        Pro_IL6 %Probability that this internal node is on.
        Pro_Ulk12 %Probability that this internal node is on.
        Pro_IL2R %Probability that this internal node is on.
        Pro_AMPK %Probability that this internal node is on.
        Pro_Tbet %Probability that this internal node is on.
        Pro_IL21 %Probability that this internal node is on.
        Pro_IL4R %Probability that this internal node is on.
        Pro_IL18R %Probability that this internal node is on.
        Pro_HPK1 %Probability that this internal node is on.
        Pro_STAT4 %Probability that this internal node is on.
        Pro_mTORC2 %Probability that this internal node is on.
        Pro_TGFbR %Probability that this internal node is on.
        Pro_ACAD %Probability that this internal node is on.
        Pro_Lip_syn %Probability that this internal node is on.
        Pro_IL4 %Probability that this internal node is on.
        Pro_AICD1 %Probability that this internal node is on.
        Pro_BCl2 %Probability that this internal node is on.
        Pro_IFNg %Probability that this internal node is on.
        Pro_STAT3 %Probability that this internal node is on.
        Pro_LXRs %Probability that this internal node is on.
        Pro_STAT1 %Probability that this internal node is on.
        Pro_STAT6 %Probability that this internal node is on.
        Pro_CFLIP %Probability that this internal node is on.
        Pro_NFkB %Probability that this internal node is on.
        Pro_Cycle %Probability that this internal node is on.
        Pro_SOCS1 %Probability that this internal node is on.
        Pro_Mem %Probability that this internal node is on.
        Pro_IRAK %Probability that this internal node is on.
        Pro_GATA3 %Probability that this internal node is on.
        Pro_IL17 %Probability that this internal node is on.
        Pro_AICD2 %Probability that this internal node is on.
        Pro_aatrans %Probability that this internal node is on.
        Pro_IL12R %Probability that this internal node is on.
        Pro_IFNgR %Probability that this internal node is on.
        Pro_RORgt %Probability that this internal node is on.
        Pro_IL2 %Probability that this internal node is on.
        Pro_aGly %Probability that this internal node is on.
        Pro_SREBP2 %Probability that this internal node is on.
        Pro_NFAT %Probability that this internal node is on.
        Pro_PI3K %Probability that this internal node is on.
        Pro_p73 %Probability that this internal node is on.
        Pro_Foxp3 %Probability that this internal node is on.
        Pro_Auto %Probability that this internal node is on.
        Pro_ERRa %Probability that this internal node is on.
        Pro_Glu_uptake %Probability that this internal node is on.
        Pro_Lip_eff %Probability that this internal node is on.
        Pro_ROS %Probability that this internal node is on.
        Pro_STAT5 %Probability that this internal node is on.
        Pro_mTORC1 %Probability that this internal node is on.
        Pro_Mito_ox %Probability that this internal node is on.
        Pro_IL21R %Probability that this internal node is on.
        BM_DNA %Amount of biomass minus DNA (g).
        DNA %Amount of DNA (g).
        AMP %Concentration of AMP (M).
        ATP %Concentration of ATP (M).                
        G1_BM_DNA %Production rate of biomass minus DNA (g per hour).
        G1_DNA %Production rate of DNA (g per hour).
        G1_ATP %Production rate of ATP (M per hour).
        G1_AMP %Production rate of AMP (M per hour).
        S_BM_DNA %Production rate of biomass minus DNA (g per hour).
        S_DNA %Production rate of DNA (g per hour).
        S_ATP %Production rate of ATP (M per hour).
        S_AMP %Production rate of AMP (M per hour).
        
    end
    
    methods
        
        function obj=CD4(ID, Location, step_size, k_tot, k_toln, k_toc)
        %This function creates an instance of the generic CD4+ T cell (an
        %object). The instance is a naive CD4+ T cell. Some of its
        %properties like its location are provided by the user. The others
        %take their default values.
            BM_DNA_initial=(64.5-8)*1e-12; %The Measurement of Lymphocyte Volume: Importance of Reference Particle Deformability and Counting Solution Tonicity (Segel et. al. 1981). Also, according to The Molecular Composition of Cells (2nd edition), water accounts for 70 % or more of total cell mass.
            DNA_initial=8e-12; %Nuclear DNA Content Varies with Cell Size across Human Cell Types (Gillooly et. al. 2015).
            AMP_initial=85e-9; %AMP-Activated Protein Kinase: the Energy Charge Hypothesis Revisited (Hardie 2001).
            ATP_initial=8500e-9; %Monitoring of ATP Levels in Red Blood Cells and T Cells of Healthy and Ill Subjects and the Effects of Age on Mitochondrial Potential (Mikirova et. al. 2004).
            div_lim_value=10; %Number of times an effector cell can divide. An Overview of the Immune System (Parkin 2001).
            mig_thres_value=6; %An effector cell can exit LN only after six or more divisions. CD4 Effector T Cell Subsets in the Response to Influenza (Roman 2002).
            act_thres_naive=24; %Number of TCR stimuli required to activate a naive cell. Two Distinct Stages in the Transition from Naive CD4 T Cells to Effectors, Early Antigen-dependent and Late Cytokine-driven Expansion and Differentiation (Jelley-Gibbs 2000).
            if nargin>0
                obj.ID=ID;
                obj.step_size=step_size;
                obj.div_lim=div_lim_value;
                obj.mig_thres=mig_thres_value;
                obj.act_thres=act_thres_naive/step_size;
                obj.State=1;
                obj.BN_count=0; %Clearly, the BN model has not been run for the newly created cell yet.
                obj.Phenotype=0;
                obj.Phase=1; %Upon activation, the CD4+ T cell enters the G1 phase where it produces biomass minus DNA.
                obj.Location=Location;
                obj.k_tot=k_tot;
                obj.k_toln=k_toln;
                obj.k_toc=k_toc;                
                obj.Divisions=0;                
                obj.Plasticity=1; %At the beginning, the cell is plastic and can transdifferentiate. Functional and Phenotypic Plasticity of CD4+ T Cell Subsets (Caza and Landas 2015).
                obj.TCR_Stimuli=0;
                obj.TCR_Strength=0;
                obj.CD28_Strength=0;
                obj.mTORC1_t=0; %At the beginning, mTORC1 is inactive.
                obj.Ribosome=0; %At the beginning, ribosome is inactive.
                obj.Resist=1; %At the beginning, the cell is resistant to AICD. Concepts of Activated T Cell Death (Brenner 2008).
                obj.IL2_sum=0;
                obj.IL4_sum=0;
                obj.IL6_sum=0;
                obj.IL12_sum=0;
                obj.IL18_sum=0;
                obj.IL23_sum=0;
                obj.IL27_sum=0;
                obj.IFNg_sum=0;
                obj.TGFb_sum=0;
                obj.Pro_IL6R=0;
                obj.Pro_Jak1=0;
                obj.Pro_Gluta=0;
                obj.Pro_IL23R=0;
                obj.Pro_Akt=0;
                obj.Pro_cMyc=0; 
                obj.Pro_FasLR=0;
                obj.Pro_IL6=0;
                obj.Pro_Ulk12=0;
                obj.Pro_IL2R=0;
                obj.Pro_AMPK=0;
                obj.Pro_Tbet=0;
                obj.Pro_IL21=0;
                obj.Pro_IL4R=0;
                obj.Pro_IL18R=0;
                obj.Pro_HPK1=0;
                obj.Pro_STAT4=0;
                obj.Pro_mTORC2=0;
                obj.Pro_TGFbR=0;
                obj.Pro_ACAD=0;
                obj.Pro_Lip_syn=0;
                obj.Pro_IL4=0;
                obj.Pro_AICD1=0;
                obj.Pro_BCl2=0;
                obj.Pro_IFNg=0;
                obj.Pro_STAT3=0;
                obj.Pro_LXRs=0;
                obj.Pro_STAT1=0;
                obj.Pro_STAT6=0;
                obj.Pro_CFLIP=0;
                obj.Pro_NFkB=0;
                obj.Pro_Cycle=0;
                obj.Pro_SOCS1=0;
                obj.Pro_Mem=0;
                obj.Pro_IRAK=0;
                obj.Pro_GATA3=0;
                obj.Pro_IL17=0;
                obj.Pro_AICD2=0;
                obj.Pro_aatrans=0;
                obj.Pro_IL12R=0;
                obj.Pro_IFNgR=0;
                obj.Pro_RORgt=0;
                obj.Pro_IL2=0;
                obj.Pro_aGly=0;
                obj.Pro_SREBP2=0;
                obj.Pro_NFAT=0;
                obj.Pro_PI3K=0;
                obj.Pro_p73=0;
                obj.Pro_Foxp3=0;
                obj.Pro_Auto=0;
                obj.Pro_ERRa=0;
                obj.Pro_Glu_uptake=0;
                obj.Pro_Lip_eff=0;
                obj.Pro_ROS=0;
                obj.Pro_STAT5=0;
                obj.Pro_mTORC1=0;
                obj.Pro_Mito_ox=0;
                obj.Pro_IL21R=0;
                obj.BM_DNA=BM_DNA_initial;
                obj.DNA=DNA_initial;
                obj.AMP=AMP_initial;
                obj.ATP=ATP_initial;
                obj.G1_BM_DNA=0; %Assume that a naive CD4+ T cell is metabolically at rest.
                obj.G1_DNA=0; %Assume that a naive CD4+ T cell is metabolically at rest.
                obj.G1_ATP=0; %Assume that a naive CD4+ T cell is metabolically at rest.
                obj.G1_AMP=0; %Assume that a naive CD4+ T cell is metabolically at rest.
                obj.S_BM_DNA=0; %Assume that a naive CD4+ T cell is metabolically at rest.
                obj.S_DNA=0; %Assume that a naive CD4+ T cell is metabolically at rest.
                obj.S_ATP=0; %Assume that a naive CD4+ T cell is metabolically at rest.
                obj.S_AMP=0; %Assume that a naive CD4+ T cell is metabolically at rest.
            end
        end
        
        function BN=BN(CD4, C_TCR, C_CD28, C_IL18, C_IFNg, C_TGFb, C_IL27, C_IL4, C_IL6, C_IL12, C_IL23, C_IL2)
        %This function takes a CD4+ T cell and various inputs and runs the
        %BN model in order to update the internal state of the cell. The
        %decision to override the initial internal state depends on the
        %plasticity of the cell. Then, assign a phenotype and update the nu
        %mber of times mTORC1 has been on so far.
            %First, convert the cytokine concentrations into the activity levels of the external nodes in the BN model.
            %Assume that the diameter of a CD4+ T cell is 10 microns to calculate its volume and hence the number of each cytokine in that volume.
            %The numbers of cytokine receptors on the cell's surface range from hundreds to thousands, 100 to 6000. Cytokine Receptors: Structure and Signal Transduction (Foxwell et. al. 1992) and Cytokine Receptor Expression on Hematopoietic Stem and Progenitor Cells (McKinstry et. al. 1997).
            %At saturation, there are around 1500 IL2/IL2R complexes on a T cell's surface. An increase in the number of complexes cannot trigger a higher proliferation rate. Computational Model for Effects of Ligand/Receptor Binding Properties on IL2 Trafficking Dynamics and T Cell Proliferation Response (Fallon et. al. 2000).        
            cell_volume=4/3*pi*(5e-5)^3;
            rep_min=100;
            rep_max=6000;
            sat_count=1500;
            N_IL18=C_IL18*cell_volume*6.02214e23;
            R_IL18=rep_min+(rep_max-rep_min)*rand;
            C_IL18=min(N_IL18, R_IL18)/sat_count;
            N_IFNg=C_IFNg*cell_volume*6.02214e23;
            R_IFNg=rep_min+(rep_max-rep_min)*rand;
            C_IFNg=min(N_IFNg, R_IFNg)/sat_count; 
            N_TGFb=C_TGFb*cell_volume*6.02214e23;
            R_TGFb=rep_min+(rep_max-rep_min)*rand;
            C_TGFb=min(N_TGFb, R_TGFb)/sat_count;
            N_IL27=C_IL27*cell_volume*6.02214e23;
            R_IL27=rep_min+(rep_max-rep_min)*rand;
            C_IL27=min(N_IL27, R_IL27)/sat_count;    
            N_IL4=C_IL4*cell_volume*6.02214e23;
            R_IL4=rep_min+(rep_max-rep_min)*rand;
            C_IL4=min(N_IL4, R_IL4)/sat_count; 
            N_IL6=C_IL6*cell_volume*6.02214e23;
            R_IL6=rep_min+(rep_max-rep_min)*rand;
            C_IL6=min(N_IL6, R_IL6)/sat_count;    
            N_IL12=C_IL12*cell_volume*6.02214e23;
            R_IL12=rep_min+(rep_max-rep_min)*rand;
            C_IL12=min(N_IL12, R_IL12)/sat_count;   
            N_IL23=C_IL23*cell_volume*6.02214e23;
            R_IL23=rep_min+(rep_max-rep_min)*rand;
            C_IL23=min(N_IL23, R_IL23)/sat_count;   
            N_IL2=C_IL2*cell_volume*6.02214e23;
            R_IL2=rep_min+(rep_max-rep_min)*rand;
            C_IL2=min(N_IL2, R_IL2)/sat_count;     
            %Second, convert mTORC1_t into a BN input. The effect of mTORC1_t on the BN model is due to chronic activation of mTORC1, so it takes time to kick in.
            %Sestrin as a feedback inhibitor of TOR that prevents age-related pathologies (Lee et. al. 2010).
            %The threshold of ten BN simulations is an assumption.
            BN_count_thres=10;
            if CD4.BN_count>=BN_count_thres
                dummy_mTORC1_t=CD4.mTORC1_t/CD4.BN_count;
            else
                dummy_mTORC1_t=0;
            end
            %Third, solve the BN model after defining iterate_run and iterate_record.
            iterate_run=250;
            iterate_record=200;
            [Temp_Pro_IL6R, Temp_Pro_Jak1, Temp_Pro_Gluta, Temp_Pro_IL23R, Temp_Pro_Akt, Temp_Pro_cMyc,...
            Temp_Pro_FasLR, Temp_Pro_IL6, Temp_Pro_Ulk12, Temp_Pro_IL2R, Temp_Pro_AMPK, Temp_Pro_Tbet,...
            Temp_Pro_IL21, Temp_Pro_IL4R, Temp_Pro_IL18R, Temp_Pro_HPK1, Temp_Pro_STAT4, Temp_Pro_mTORC2,...
            Temp_Pro_TGFbR, Temp_Pro_ACAD, Temp_Pro_Lip_syn, Temp_Pro_IL4, Temp_Pro_AICD1, Temp_Pro_BCl2,...
            Temp_Pro_IFNg, Temp_Pro_STAT3, Temp_Pro_LXRs, Temp_Pro_STAT1, Temp_Pro_STAT6, Temp_Pro_CFLIP,...
            Temp_Pro_NFkB, Temp_Pro_Cycle, Temp_Pro_SOCS1, Temp_Pro_Mem, Temp_Pro_IRAK, Temp_Pro_GATA3,...
            Temp_Pro_IL17, Temp_Pro_AICD2, Temp_Pro_aatrans, Temp_Pro_IL12R, Temp_Pro_IFNgR, Temp_Pro_RORgt,...
            Temp_Pro_IL2, Temp_Pro_aGly, Temp_Pro_SREBP2, Temp_Pro_NFAT, Temp_Pro_PI3K, Temp_Pro_p73,...
            Temp_Pro_Foxp3, Temp_Pro_Auto, Temp_Pro_ERRa, Temp_Pro_Glu_uptake, Temp_Pro_Lip_eff,...
            Temp_Pro_ROS, Temp_Pro_STAT5, Temp_Pro_mTORC1, Temp_Pro_Mito_ox, Temp_Pro_IL21R]=BN_model(C_IL18,...
            CD4.Resist, dummy_mTORC1_t, C_IFNg, C_TGFb, C_IL27, C_IL4, CD4.AMP/CD4.ATP, C_CD28,...
            C_TCR, C_IL6, CD4.Ribosome, C_IL12, C_IL23, C_IL2, CD4.Pro_IL6R, CD4.Pro_Jak1,...
            CD4.Pro_Gluta, CD4.Pro_IL23R, CD4.Pro_Akt, CD4.Pro_cMyc, CD4.Pro_FasLR, CD4.Pro_IL6,...
            CD4.Pro_Ulk12, CD4.Pro_IL2R, CD4.Pro_AMPK, CD4.Pro_Tbet, CD4.Pro_IL21, CD4.Pro_IL4R,...
            CD4.Pro_IL18R, CD4.Pro_HPK1, CD4.Pro_STAT4, CD4.Pro_mTORC2, CD4.Pro_TGFbR, CD4.Pro_ACAD,...
            CD4.Pro_Lip_syn, CD4.Pro_IL4, CD4.Pro_AICD1, CD4.Pro_BCl2, CD4.Pro_IFNg, CD4.Pro_STAT3,...
            CD4.Pro_LXRs, CD4.Pro_STAT1, CD4.Pro_STAT6, CD4.Pro_CFLIP, CD4.Pro_NFkB, CD4.Pro_Cycle,...
            CD4.Pro_SOCS1, CD4.Pro_Mem, CD4.Pro_IRAK, CD4.Pro_GATA3, CD4.Pro_IL17, CD4.Pro_AICD2,...
            CD4.Pro_aatrans, CD4.Pro_IL12R, CD4.Pro_IFNgR, CD4.Pro_RORgt, CD4.Pro_IL2, CD4.Pro_aGly,...
            CD4.Pro_SREBP2, CD4.Pro_NFAT, CD4.Pro_PI3K, CD4.Pro_p73, CD4.Pro_Foxp3, CD4.Pro_Auto,...
            CD4.Pro_ERRa, CD4.Pro_Glu_uptake, CD4.Pro_Lip_eff, CD4.Pro_ROS, CD4.Pro_STAT5, CD4.Pro_mTORC1,...
            CD4.Pro_Mito_ox, CD4.Pro_IL21R, iterate_run, iterate_record);
            %Increment the number of times the BN model has been run for this cell.
            CD4.BN_count=CD4.BN_count+1;
            %Update the cell's internal state with the output from the BN
            %model. These nodes only update when the cell is plastic.
            if rand<CD4.Plasticity
                CD4.Pro_Tbet=Temp_Pro_Tbet;
                CD4.Pro_GATA3=Temp_Pro_GATA3;
                CD4.Pro_RORgt=Temp_Pro_RORgt;
                CD4.Pro_Foxp3=Temp_Pro_Foxp3;
                CD4.Pro_IL2=Temp_Pro_IL2;
                CD4.Pro_IL4=Temp_Pro_IL4;               
                CD4.Pro_IL6=Temp_Pro_IL6;
                CD4.Pro_IL17=Temp_Pro_IL17;          
                CD4.Pro_IL21=Temp_Pro_IL21;
                CD4.Pro_IFNg=Temp_Pro_IFNg;
            end
            %Update the cell's internal state with the output from the BN
            %model. These nodes update even for terminally differentiated
            %cells. Apart from mTORC1, convert the activity levels of the
            %metabolic and attribute nodes before updating the cell. 
            %These nodes are end nodes, meaning they do not affect the
            %states of any other nodes in the BN model. When the BN model 
            %is run, their initial states will be overwritten immediately.
            CD4.Pro_IL6R=Temp_Pro_IL6R;
            CD4.Pro_Jak1=Temp_Pro_Jak1;
            CD4.Pro_IL23R=Temp_Pro_IL23R;
            CD4.Pro_Akt=Temp_Pro_Akt;
            CD4.Pro_cMyc=Temp_Pro_cMyc;
            CD4.Pro_FasLR=Temp_Pro_FasLR;
            CD4.Pro_Ulk12=Temp_Pro_Ulk12;
            CD4.Pro_IL2R=Temp_Pro_IL2R;
            CD4.Pro_AMPK=Temp_Pro_AMPK;
            CD4.Pro_IL4R=Temp_Pro_IL4R;
            CD4.Pro_IL18R=Temp_Pro_IL18R;
            CD4.Pro_HPK1=Temp_Pro_HPK1;
            CD4.Pro_STAT4=Temp_Pro_STAT4;
            CD4.Pro_mTORC2=Temp_Pro_mTORC2;
            CD4.Pro_TGFbR=Temp_Pro_TGFbR;
            CD4.Pro_BCl2=Temp_Pro_BCl2;
            CD4.Pro_STAT3=Temp_Pro_STAT3;
            CD4.Pro_LXRs=Temp_Pro_LXRs;
            CD4.Pro_STAT1=Temp_Pro_STAT1;
            CD4.Pro_STAT6=Temp_Pro_STAT6;
            CD4.Pro_CFLIP=Temp_Pro_CFLIP;
            CD4.Pro_NFkB=Temp_Pro_NFkB;            
            CD4.Pro_SOCS1=Temp_Pro_SOCS1;
            CD4.Pro_IRAK=Temp_Pro_IRAK;
            CD4.Pro_IL12R=Temp_Pro_IL12R;
            CD4.Pro_IFNgR=Temp_Pro_IFNgR;
            CD4.Pro_SREBP2=Temp_Pro_SREBP2;
            CD4.Pro_NFAT=Temp_Pro_NFAT;
            CD4.Pro_PI3K=Temp_Pro_PI3K;
            CD4.Pro_p73=Temp_Pro_p73;
            CD4.Pro_ERRa=Temp_Pro_ERRa;
            CD4.Pro_ROS=Temp_Pro_ROS;
            CD4.Pro_STAT5=Temp_Pro_STAT5;
            CD4.Pro_IL21R=Temp_Pro_IL21R;
            CD4.Pro_Gluta=Metabolism_conversion(Temp_Pro_Gluta, 6);
            CD4.Pro_Lip_syn=Metabolism_conversion(Temp_Pro_Lip_syn, 7);            
            CD4.Pro_aatrans=Metabolism_conversion(Temp_Pro_aatrans, 3);            
            CD4.Pro_aGly=Metabolism_conversion(Temp_Pro_aGly, 1);            
            CD4.Pro_Glu_uptake=Metabolism_conversion(Temp_Pro_Glu_uptake, 2);
            CD4.Pro_Lip_eff=Metabolism_conversion(Temp_Pro_Lip_eff, 5);            
            CD4.Pro_Mito_ox=Metabolism_conversion(Temp_Pro_Mito_ox, 4);
            CD4.Pro_Auto=Attribute_conversion(Temp_Pro_Auto, 2);
            CD4.Pro_AICD1=Attribute_conversion(Temp_Pro_AICD1, 5);
            CD4.Pro_AICD2=Attribute_conversion(Temp_Pro_AICD2, 6);
            CD4.Pro_ACAD=Attribute_conversion(Temp_Pro_ACAD, 4);
            CD4.Pro_Mem=Attribute_conversion(Temp_Pro_Mem, 3);
            CD4.Pro_Cycle=Attribute_conversion(Temp_Pro_Cycle, 1);
            CD4.Pro_mTORC1=Temp_Pro_mTORC1;
            %Assign a phenotype to the cell. Remember to convert the
            %activity levels of the transcription factors into their
            %empirical cumulative distribution function values. However,
            %the activity levels recorded inside the cell will not change.
            dummy_Tbet=rand;
            dummy_GATA3=rand;
            dummy_RORgt=rand;
            dummy_Foxp3=rand;
            scaled_Tbet=TF_conversion(CD4.Pro_Tbet, 1);
            scaled_GATA3=TF_conversion(CD4.Pro_GATA3, 2);
            scaled_RORgt=TF_conversion(CD4.Pro_RORgt, 3);
            scaled_Foxp3=TF_conversion(CD4.Pro_Foxp3, 4);
            if ((dummy_Tbet<scaled_Tbet) & (dummy_GATA3<scaled_GATA3) & (dummy_RORgt<scaled_RORgt) & (dummy_Foxp3<scaled_Foxp3))
                CD4.Phenotype=1;
            elseif ((dummy_GATA3<scaled_GATA3) & (dummy_RORgt<scaled_RORgt) & (dummy_Foxp3<scaled_Foxp3))
                CD4.Phenotype=2;
            elseif ((dummy_Tbet<scaled_Tbet) & (dummy_RORgt<scaled_RORgt) & (dummy_Foxp3<scaled_Foxp3))
                CD4.Phenotype=3;
            elseif ((dummy_Tbet<scaled_Tbet) & (dummy_GATA3<scaled_GATA3) & (dummy_Foxp3<scaled_Foxp3))
                CD4.Phenotype=4;
            elseif ((dummy_Tbet<scaled_Tbet) & (dummy_GATA3<scaled_GATA3) & (dummy_RORgt<scaled_RORgt))
                CD4.Phenotype=5;
            elseif ((dummy_RORgt<scaled_RORgt) & (dummy_Foxp3<scaled_Foxp3))
                CD4.Phenotype=6;
            elseif ((dummy_Tbet<scaled_Tbet) & (dummy_GATA3<scaled_GATA3))
                CD4.Phenotype=7;
            elseif ((dummy_Tbet<scaled_Tbet) & (dummy_RORgt<scaled_RORgt))
                CD4.Phenotype=8;
            elseif ((dummy_Tbet<scaled_Tbet) & (dummy_Foxp3<scaled_Foxp3))
                CD4.Phenotype=9;
            elseif ((dummy_GATA3<scaled_GATA3) & (dummy_Foxp3<scaled_Foxp3))
                CD4.Phenotype=10;
            elseif ((dummy_GATA3<scaled_GATA3) & (dummy_RORgt<scaled_RORgt))
                CD4.Phenotype=11;
            elseif (dummy_Tbet<scaled_Tbet)
                CD4.Phenotype=12; %Th1.
            elseif (dummy_GATA3<scaled_GATA3)
                CD4.Phenotype=13; %Th2.
            elseif (dummy_RORgt<scaled_RORgt)
                CD4.Phenotype=14; %Th17.
            elseif (dummy_Foxp3<scaled_Foxp3)
                CD4.Phenotype=15; %Treg.
            else
                CD4.Phenotype=0; %Th0.
            end
            %Update the number of times mTORC1 has been on so far.
            %Remember to convert the activity level of mTORC1 into its
            %empirical cumulative distribution function value. However, the
            %activity level recorded inside the cell will not change.
            scaled_mTORC1=Attribute_conversion(CD4.Pro_mTORC1, 7);
            if rand<scaled_mTORC1
                CD4.mTORC1_t=CD4.mTORC1_t+1;
            end
            %Return the updated cell.
            BN=CD4;
        end
        
        function [MN, MN_input, MN_output]=MN(CD4, model, model_phenotype, protein_max, solution1, solution2, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, MN_input, MN_output)
        %This function takes a CD4+ T cell and solves a suitable MN
        %model under the constraints of the internal state of the cell.
        %Then, it extracts various production rates from the solution of
        %the MN model. Each of MN_output1 and MN_output2 contains the
        %production rates of biomass, biomass minus DNA, DNA, protein, ATP,
        %and AMP. Because linear programming is an expensive step, the
        %function will search for suitable MN_output1 and MN_output2 values
        %in the library of pre-calculated values. If the library has
        %nothing suitable for the internal state of the cell, the function
        %will perform the linear programming operation and expand the
        %library after that.
            %First, reset the variables.
            dummy=size(MN_input, 1);
            best_index=0;
            best_distance=0;
            entry_range=0.05;
            %Consider each entry in the library.
            for i=1:dummy
                score=0;
                %First, check if the activity levels of the metabolic nodes
                %in the current library entry are sufficiently close to the
                %internal state of the cell.
                if ((CD4.Pro_aGly-entry_range)<MN_input(i, 1)) & (MN_input(i, 1)<(CD4.Pro_aGly+entry_range))
                    score=score+1;
                end
                if ((CD4.Pro_Glu_uptake-entry_range)<MN_input(i, 2)) & (MN_input(i, 2)<(CD4.Pro_Glu_uptake+entry_range))
                    score=score+1;
                end
                if ((CD4.Pro_aatrans-entry_range)<MN_input(i, 3)) & (MN_input(i, 3)<(CD4.Pro_aatrans+entry_range))
                    score=score+1;
                end
                if ((CD4.Pro_Mito_ox-entry_range)<MN_input(i, 4)) & (MN_input(i, 4)<(CD4.Pro_Mito_ox+entry_range))
                    score=score+1;
                end
                if ((CD4.Pro_Lip_eff-entry_range)<MN_input(i, 5)) & (MN_input(i, 5)<(CD4.Pro_Lip_eff+entry_range))
                    score=score+1;
                end
                if ((CD4.Pro_Gluta-entry_range)<MN_input(i, 6)) & (MN_input(i, 6)<(CD4.Pro_Gluta+entry_range))
                    score=score+1;
                end
                if ((CD4.Pro_Lip_syn-entry_range)<MN_input(i, 7)) & (MN_input(i, 7)<(CD4.Pro_Lip_syn+entry_range))
                    score=score+1;
                end
                %If the entry is sufficiently close to the internal state
                %of the cell in all 7 dimensions, calculate the overall
                %distance.
                if score==7
                    %Each term is the relative difference between the entry
                    %value of a metabolic node and the internal state.
                    distance=abs(MN_input(i, 1)-CD4.Pro_aGly)/CD4.Pro_aGly+abs(MN_input(i, 2)-CD4.Pro_Glu_uptake)/CD4.Pro_Glu_uptake+abs(MN_input(i, 3)-CD4.Pro_aatrans)/CD4.Pro_aatrans+abs(MN_input(i, 4)-CD4.Pro_Mito_ox)/CD4.Pro_Mito_ox+abs(MN_input(i, 5)-CD4.Pro_Lip_eff)/CD4.Pro_Lip_eff+abs(MN_input(i, 6)-CD4.Pro_Gluta)/CD4.Pro_Gluta+abs(MN_input(i, 7)-CD4.Pro_Lip_syn)/CD4.Pro_Lip_syn;
                    %Check if the distance is the smallest calculated so
                    %far.
                    if best_index==0|distance<best_distance
                        best_distance=distance;
                        best_index=i;
                    end
                end
            end
            %%Assume that the diameter of a CD4+ T cell is 10 microns to calculate its volume.
            cell_volume=4/3*pi*(5e-5)^3;
            %After scanning through the whole library, check if there are
            %any entries that are sufficiently close to the internal state
            %of the cell.
            if best_index~=0
                %There is at least one entry that is sufficiently close to
                %the internal state of the cell. Pick the entry with the
                %smallest distance and use that as MN_output1 and
                %MN_output2.
                MN_output1=[MN_output(best_index, 1), MN_output(best_index, 2), MN_output(best_index, 3), MN_output(best_index, 4), MN_output(best_index, 5), MN_output(best_index, 6)];
                MN_output2=[MN_output(best_index, 7), MN_output(best_index, 8), MN_output(best_index, 9), MN_output(best_index, 10), MN_output(best_index, 11), MN_output(best_index, 12)];
                %In G1, the cell optimises the production of biomass
                %minus DNA, so the production rate of DNA is zero.
                %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                Library_G1_BM_DNA=(MN_output1(2)/log(2))*CD4.BM_DNA;
                Library_G1_DNA=0;
                Library_G1_ATP=MN_output1(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                Library_G1_AMP=MN_output1(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                %In S, the cell optimises the production of DNA.
                %However, we cannot convert the molar production rate
                %of DNA to a mass production rate, so we split the
                %biomass production rate into two parts.
                %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                Library_S_BM_DNA=(MN_output2(1)/log(2))*CD4.BM_DNA;
                Library_S_DNA=(MN_output2(1)/log(2))*CD4.DNA;
                Library_S_ATP=MN_output2(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                Library_S_AMP=MN_output2(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                %Average the optimised protein production rates and use the
                %result to get the activity level of ribosome.
                Library_Ribosome=(MN_output1(4)+MN_output2(4))/protein_max/2;
                %Return the production rates.
                MN=[Library_G1_BM_DNA, Library_G1_DNA, Library_G1_ATP, Library_G1_AMP, Library_S_BM_DNA, Library_S_DNA, Library_S_ATP, Library_S_AMP, Library_Ribosome];                
            else
                %Since there is not a single entry in the library that is
                %close enough to the internal state, it is necessary to
                %calculate directly.
                if model_phenotype==0 %Th0.
                    %Optimise the MN model. Optimisation must be carried out twice, once for phase 1 (G1) and once for phase 2 (S).
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    %Some of the BN outputs are normalised before being fed into the metabolic models because they are not used as probabilities.
                    %In other words, there is a chance that aerobic glycolysis is maximum when the controlling node is on 25 % of the time only.
                    MN_output1 = Th0_MN_control(model, 1, solution1, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    MN_output2 = Th0_MN_control(model, 2, solution2, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    %Record the new entry in the library.
                    MN_input_new=[CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn];
                    MN_input=[MN_input; MN_input_new];
                    MN_output_new=[MN_output1, MN_output2];
                    MN_output=[MN_output; MN_output_new];
                    %In G1, the cell optimises the production of biomass
                    %minus DNA, so the production rate of DNA is zero.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Th0_G1_BM_DNA=(MN_output1(2)/log(2))*CD4.BM_DNA;
                    Th0_G1_DNA=0;
                    Th0_G1_ATP=MN_output1(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Th0_G1_AMP=MN_output1(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %In S, the cell optimises the production of DNA.
                    %However, we cannot convert the molar production rate
                    %of DNA to a mass production rate, so we split the
                    %biomass production rate into two parts.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Th0_S_BM_DNA=(MN_output2(1)/log(2))*CD4.BM_DNA;
                    Th0_S_DNA=(MN_output2(1)/log(2))*CD4.DNA;
                    Th0_S_ATP=MN_output2(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Th0_S_AMP=MN_output2(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %Average the optimised protein production rates and use the
                    %result to get the activity level of ribosome.
                    Th0_Ribosome=(MN_output1(4)+MN_output2(4))/protein_max/2;
                    %Return the production rates for Th0.
                    MN=[Th0_G1_BM_DNA, Th0_G1_DNA, Th0_G1_ATP, Th0_G1_AMP, Th0_S_BM_DNA, Th0_S_DNA, Th0_S_ATP, Th0_S_AMP, Th0_Ribosome];
                elseif model_phenotype==12 %Th1.
                    %Optimise the MN model. Optimisation must be carried out twice, once for phase 1 (G1) and once for phase 2 (S).
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    %Some of the BN outputs are normalised before being fed into the metabolic models because they are not used as probabilities.
                    %In other words, there is a chance that aerobic glycolysis is maximum when the controlling node is on 25 % of the time only.                
                    MN_output1 = Th1_MN_control(model, 1, solution1, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    MN_output2 = Th1_MN_control(model, 2, solution2, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    %Record the new entry in the library.                    
                    MN_input_new=[CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn];
                    MN_input=[MN_input; MN_input_new];
                    MN_output_new=[MN_output1, MN_output2];
                    MN_output=[MN_output; MN_output_new];                    
                    %In G1, the cell optimises the production of biomass
                    %minus DNA, so the production rate of DNA is zero.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Th1_G1_BM_DNA=(MN_output1(2)/log(2))*CD4.BM_DNA;
                    Th1_G1_DNA=0;
                    Th1_G1_ATP=MN_output1(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Th1_G1_AMP=MN_output1(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %In S, the cell optimises the production of DNA.
                    %However, we cannot convert the molar production rate
                    %of DNA to a mass production rate, so we split the
                    %biomass production rate into two parts.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Th1_S_BM_DNA=(MN_output2(1)/log(2))*CD4.BM_DNA;
                    Th1_S_DNA=(MN_output2(1)/log(2))*CD4.DNA;
                    Th1_S_ATP=MN_output2(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Th1_S_AMP=MN_output2(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %Average the optimised protein production rates and use the
                    %result to get the activity level of ribosome.
                    Th1_Ribosome=(MN_output1(4)+MN_output2(4))/protein_max/2;
                    %Return the production rates for Th1.
                    MN=[Th1_G1_BM_DNA, Th1_G1_DNA, Th1_G1_ATP, Th1_G1_AMP, Th1_S_BM_DNA, Th1_S_DNA, Th1_S_ATP, Th1_S_AMP, Th1_Ribosome];             
                elseif model_phenotype==13 %Th2.  
                    %Optimise the MN model. Optimisation must be carried out twice, once for phase 1 (G1) and once for phase 2 (S).
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    %Some of the BN outputs are normalised before being fed into the metabolic models because they are not used as probabilities.
                    %In other words, there is a chance that aerobic glycolysis is maximum when the controlling node is on 25 % of the time only.
                    MN_output1 = Th2_MN_control(model, 1, solution1, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    MN_output2 = Th2_MN_control(model, 2, solution2, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    %Record the new entry in the library.                    
                    MN_input_new=[CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn];
                    MN_input=[MN_input; MN_input_new];
                    MN_output_new=[MN_output1, MN_output2];
                    MN_output=[MN_output; MN_output_new];                    
                    %In G1, the cell optimises the production of biomass
                    %minus DNA, so the production rate of DNA is zero.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Th2_G1_BM_DNA=(MN_output1(2)/log(2))*CD4.BM_DNA;
                    Th2_G1_DNA=0;
                    Th2_G1_ATP=MN_output1(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Th2_G1_AMP=MN_output1(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %In S, the cell optimises the production of DNA.
                    %However, we cannot convert the molar production rate
                    %of DNA to a mass production rate, so we split the
                    %biomass production rate into two parts.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Th2_S_BM_DNA=(MN_output2(1)/log(2))*CD4.BM_DNA;
                    Th2_S_DNA=(MN_output2(1)/log(2))*CD4.DNA;
                    Th2_S_ATP=MN_output2(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Th2_S_AMP=MN_output2(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %Average the optimised protein production rates and use the
                    %result to get the activity level of ribosome.
                    Th2_Ribosome=(MN_output1(4)+MN_output2(4))/protein_max/2;
                    %Return the production rates for Th2.
                    MN=[Th2_G1_BM_DNA, Th2_G1_DNA, Th2_G1_ATP, Th2_G1_AMP, Th2_S_BM_DNA, Th2_S_DNA, Th2_S_ATP, Th2_S_AMP, Th2_Ribosome];              
                elseif model_phenotype==14 %Th17. 
                    %Optimise the MN model. Optimisation must be carried out twice, once for phase 1 (G1) and once for phase 2 (S).
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    %Some of the BN outputs are normalised before being fed into the metabolic models because they are not used as probabilities.
                    %In other words, there is a chance that aerobic glycolysis is maximum when the controlling node is on 25 % of the time only.                
                    MN_output1 = Th17_MN_control(model, 1, solution1, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    MN_output2 = Th17_MN_control(model, 2, solution2, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    %Record the new entry in the library.                    
                    MN_input_new=[CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn];
                    MN_input=[MN_input; MN_input_new];
                    MN_output_new=[MN_output1, MN_output2];
                    MN_output=[MN_output; MN_output_new];                    
                    %In G1, the cell optimises the production of biomass
                    %minus DNA, so the production rate of DNA is zero.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Th17_G1_BM_DNA=(MN_output1(2)/log(2))*CD4.BM_DNA;
                    Th17_G1_DNA=0;
                    Th17_G1_ATP=MN_output1(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Th17_G1_AMP=MN_output1(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %In S, the cell optimises the production of DNA.
                    %However, we cannot convert the molar production rate
                    %of DNA to a mass production rate, so we split the
                    %biomass production rate into two parts.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Th17_S_BM_DNA=(MN_output2(1)/log(2))*CD4.BM_DNA;
                    Th17_S_DNA=(MN_output2(1)/log(2))*CD4.DNA;
                    Th17_S_ATP=MN_output2(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Th17_S_AMP=MN_output2(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %Average the optimised protein production rates and use the
                    %result to get the activity level of ribosome.
                    Th17_Ribosome=(MN_output1(4)+MN_output2(4))/protein_max/2;
                    %Return the production rates for Th17.
                    MN=[Th17_G1_BM_DNA, Th17_G1_DNA, Th17_G1_ATP, Th17_G1_AMP, Th17_S_BM_DNA, Th17_S_DNA, Th17_S_ATP, Th17_S_AMP, Th17_Ribosome];             
                elseif model_phenotype==15 %Treg.
                    %Optimise the MN model. Optimisation must be carried out twice, once for phase 1 (G1) and once for phase 2 (S).
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    %Some of the BN outputs are normalised before being fed into the metabolic models because they are not used as probabilities.
                    %In other words, there is a chance that aerobic glycolysis is maximum when the controlling node is on 25 % of the time only.                
                    MN_output1 = Treg_MN_control(model, 1, solution1, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    MN_output2 = Treg_MN_control(model, 2, solution2, aa_ex, aa_tran, glu_ex, glu_tran, glutaminolysis, glycolysis, lip_ex, lip_ox, lip_syn, ox_phos, pyruvate_into_mito, ATP_fluxes, ATP_coeff, AMP_fluxes, AMP_coeff, CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn);
                    %Record the new entry in the library.                    
                    MN_input_new=[CD4.Pro_aGly, CD4.Pro_Glu_uptake, CD4.Pro_aatrans, CD4.Pro_Mito_ox, CD4.Pro_Lip_eff, CD4.Pro_Gluta, CD4.Pro_Lip_syn];
                    MN_input=[MN_input; MN_input_new];
                    MN_output_new=[MN_output1, MN_output2];
                    MN_output=[MN_output; MN_output_new];                    
                    %In G1, the cell optimises the production of biomass
                    %minus DNA, so the production rate of DNA is zero.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Treg_G1_BM_DNA=(MN_output1(2)/log(2))*CD4.BM_DNA;
                    Treg_G1_DNA=0;
                    Treg_G1_ATP=MN_output1(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Treg_G1_AMP=MN_output1(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %In S, the cell optimises the production of DNA.
                    %However, we cannot convert the molar production rate
                    %of DNA to a mass production rate, so we split the
                    %biomass production rate into two parts.
                    %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
                    Treg_S_BM_DNA=(MN_output2(1)/log(2))*CD4.BM_DNA;
                    Treg_S_DNA=(MN_output2(1)/log(2))*CD4.DNA;
                    Treg_S_ATP=MN_output2(5)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    Treg_S_AMP=MN_output2(6)*(CD4.BM_DNA+CD4.DNA)/1000/cell_volume;
                    %Average the optimised protein production rates and use the
                    %result to get the activity level of ribosome.
                    Treg_Ribosome=(MN_output1(4)+MN_output2(4))/protein_max/2;
                    %Return the production rates for Treg and the updated library.
                    MN=[Treg_G1_BM_DNA, Treg_G1_DNA, Treg_G1_ATP, Treg_G1_AMP, Treg_S_BM_DNA, Treg_S_DNA, Treg_S_ATP, Treg_S_AMP, Treg_Ribosome];
                end
            end
        end
        
        function metab=metab(CD4)
        %This function takes a cell as the input, checks its phase, and
        %updates its content accordingly. Then, the function decides whether it undergoes
        %autophagy. Finally, based on the updated content, it updates its
        %phase. Note that, at the moment, there is not a model to describe
        %autophagy, so the answer as to whether the cell undergoes
        %autophagy is always 'no'.
            %Check the cell's phase and update its content using the production rates appropriate for the phase.
            %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
            %Make sure none of the mass or concentration values are negative.
            if CD4.Phase==1
                dummy_BM_DNA=CD4.BM_DNA+CD4.G1_BM_DNA*CD4.step_size;
                if dummy_BM_DNA>0
                    CD4.BM_DNA=dummy_BM_DNA;
                else
                    CD4.BM_DNA=0;
                end
                dummy_DNA=CD4.DNA+CD4.G1_DNA*CD4.step_size;
                if dummy_DNA>0
                    CD4.DNA=dummy_DNA;
                else
                    CD4.DNA=0;
                end
                dummy_ATP=CD4.ATP+CD4.G1_ATP*CD4.step_size;
                if dummy_ATP>0
                    CD4.ATP=dummy_ATP;
                else
                    CD4.ATP=0;
                end
                dummy_AMP=CD4.AMP+CD4.G1_AMP*CD4.step_size;
                if dummy_AMP>0
                    CD4.AMP=dummy_AMP;
                else
                    CD4.AMP=0;
                end
            elseif CD4.Phase==2
                dummy_BM_DNA=CD4.BM_DNA+CD4.S_BM_DNA*CD4.step_size;
                if dummy_BM_DNA>0
                    CD4.BM_DNA=dummy_BM_DNA;
                else
                    CD4.BM_DNA=0;
                end
                dummy_DNA=CD4.DNA+CD4.S_DNA*CD4.step_size;
                if dummy_DNA>0
                    CD4.DNA=dummy_DNA;
                else
                    CD4.DNA=0;
                end
                dummy_ATP=CD4.ATP+CD4.S_ATP*CD4.step_size;
                if dummy_ATP>0
                    CD4.ATP=dummy_ATP;
                else
                    CD4.ATP=0;
                end             
                dummy_AMP=CD4.AMP+CD4.S_AMP*CD4.step_size;
                if dummy_AMP>0
                    CD4.AMP=dummy_AMP;
                else
                    CD4.AMP=0;
                end
            end
            %Define the initial values of BM_DNA, DNA, AMP, and ATP.
            BM_DNA_initial=(64.5-8)*1e-12; %The Measurement of Lymphocyte Volume: Importance of Reference Particle Deformability and Counting Solution Tonicity (Segel et. al. 1981). Also, according to The Molecular Composition of Cells (2nd edition), water accounts for 70 % or more of total cell mass.
            DNA_initial=8e-12; %Nuclear DNA Content Varies with Cell Size across Human Cell Types (Gillooly et. al. 2015).
            AMP_initial=85e-9; %AMP-Activated Protein Kinase: the Energy Charge Hypothesis Revisited (Hardie 2001).
            ATP_initial=8500e-9; %Monitoring of ATP Levels in Red Blood Cells and T Cells of Healthy and Ill Subjects and the Effects of Age on Mitochondrial Potential (Mikirova et. al. 2004).
            %If the cell undergoes autophagy, just reset its content to the
            %initial values. Note that this is a simplification because a
            %whole new intracellular model must be built to describe the
            %details of autophagy. However, this is too arbitrary, so for
            %the moment, the cell is programmed to not undergo autophagy.
            if rand<CD4.Pro_Auto*0
                CD4.BM_DNA=BM_DNA_initial;
                CD4.DNA=DNA_initial;
                CD4.AMP=AMP_initial;
                CD4.ATP=ATP_initial;
            end
            %Update the cell's phase. If the cell has enough biomass minus
            %DNA and DNA, it is ready for division, so it enters phase 3
            %(G2). If it has enough biomass minus DNA but not enough DNA,
            %it enters phase 2 (S). If it needs more of both or more of 
            %biomass minus DNA, it stays in phase 1 (G1).
            %The Mammalian Cell Cycle: an Overview (Harper and Brooks 2005).
            if CD4.BM_DNA >= 2*BM_DNA_initial & CD4.DNA >= 2*DNA_initial
                CD4.Phase=3;
            elseif CD4.BM_DNA >= 2*BM_DNA_initial
                CD4.Phase=2;
            else
                CD4.Phase=1;
            end
            %Return the updated cell.
            metab=CD4;
        end

    end
    
end