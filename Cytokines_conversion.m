function AL_scaled=Cytokines_conversion(AL, label)
%This function takes the activity level of a cytokine production node and normalises it.
    if label==1 %IL2.
        AL_scaled=AL/0.925;
    elseif label==2 %IL4.
        AL_scaled=AL/0.84;
    elseif label==3 %IL6.
        AL_scaled=AL;
    elseif label==4 %IL17.
        AL_scaled=AL;
    elseif label==5 %IL21.
        AL_scaled=AL;
    elseif label==6 %IFNg.
        AL_scaled=AL;     
    end
end