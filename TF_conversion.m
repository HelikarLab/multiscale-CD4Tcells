function AL_scaled=TF_conversion(AL, label)
%This function takes the activity level of a transcription factor node and normalises it.
    if label==1 %Tbet.
        AL_scaled=AL;
    elseif label==2 %GATA3.
        AL_scaled=AL;
    elseif label==3 %RORgt.
        AL_scaled=AL;
    elseif label==4 %Foxp3.
        AL_scaled=AL/0.51;
    end
end