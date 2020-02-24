function AL_scaled=Attribute_conversion(AL, label)
%This function takes the activity level of an attribute node and normalises it.
    if label==1 %Cycle.
        AL_scaled=AL;
    elseif label==2 %Auto.
        AL_scaled=AL/0.5;
    elseif label==3 %Mem.
        AL_scaled=AL/0.975;
    elseif label==4 %ACAD.
        AL_scaled=0.025+AL;
    elseif label==5 %AICD1.
        AL_scaled=AL/0.975;
    elseif label==6 %AICD2.
        AL_scaled=AL;
    elseif label==7 %mTORC1.
        AL_scaled=AL/0.97;
    end
end