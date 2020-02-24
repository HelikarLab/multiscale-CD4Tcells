function AL_scaled=Metabolism_conversion(AL, label)
%This function takes the activity level of a metabolic node and if it is
%within its random sample (minus the zeros), estimates the value given by
%the empirical cumulative distribution function, and uses it as the converted
%activity level. If it is smaller than or larger than the random sample
%range, use the original activity level or one respectively. Note that the
%empirical cumulative distribution functions are based on the random
%samples WITHOUT the zeros.
    if label==1 %aGly.
        if AL<0.005
            AL_scaled=AL;
        elseif AL>0.925
            AL_scaled=1;
        else
            coeff=[0, 0, 0, 0, -121.62, 373.07, -448.98, 268.32, -83.046, 12.569, 0.2524];
            AL_scaled=coeff(1)*(AL^10)+coeff(2)*(AL^9)+coeff(3)*(AL^8)+coeff(4)*(AL^7)+coeff(5)*(AL^6)+coeff(6)*(AL^5)+coeff(7)*(AL^4)+coeff(8)*(AL^3)+coeff(9)*(AL^2)+coeff(10)*(AL)+coeff(11);
        end
    elseif label==2 %Glu_uptake.
        if AL<0.005
            AL_scaled=AL;
        else
            coeff=[0, 0, 0, 0, -40.9, 138.09, -183.9, 123.22, -44.353, 8.7181, 0.11538];
            AL_scaled=coeff(1)*(AL^10)+coeff(2)*(AL^9)+coeff(3)*(AL^8)+coeff(4)*(AL^7)+coeff(5)*(AL^6)+coeff(6)*(AL^5)+coeff(7)*(AL^4)+coeff(8)*(AL^3)+coeff(9)*(AL^2)+coeff(10)*(AL)+coeff(11);
        end
    elseif label==3 %aatrans.
        if AL<0.005
            AL_scaled=AL;
        else
            coeff=[0, 0, 0, 0, -30.53, 103.9, -140.1, 95.927, -36.091, 7.8139, 0.066028];
            AL_scaled=coeff(1)*(AL^10)+coeff(2)*(AL^9)+coeff(3)*(AL^8)+coeff(4)*(AL^7)+coeff(5)*(AL^6)+coeff(6)*(AL^5)+coeff(7)*(AL^4)+coeff(8)*(AL^3)+coeff(9)*(AL^2)+coeff(10)*(AL)+coeff(11);
        end
    elseif label==4 %Mito_ox.
        if AL<0.005
            AL_scaled=AL;
        else
            coeff=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0.99991, 0.00012598];
            AL_scaled=coeff(1)*(AL^10)+coeff(2)*(AL^9)+coeff(3)*(AL^8)+coeff(4)*(AL^7)+coeff(5)*(AL^6)+coeff(6)*(AL^5)+coeff(7)*(AL^4)+coeff(8)*(AL^3)+coeff(9)*(AL^2)+coeff(10)*(AL)+coeff(11);
        end
    elseif label==5 %Lip_eff.
        if AL<0.005
            AL_scaled=AL;
        else
            coeff=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0.99994, -0.0001152];
            AL_scaled=coeff(1)*(AL^10)+coeff(2)*(AL^9)+coeff(3)*(AL^8)+coeff(4)*(AL^7)+coeff(5)*(AL^6)+coeff(6)*(AL^5)+coeff(7)*(AL^4)+coeff(8)*(AL^3)+coeff(9)*(AL^2)+coeff(10)*(AL)+coeff(11);
        end        
    elseif label==6 %Gluta.
        if AL<0.005
            AL_scaled=AL;
        else
            coeff=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0.99991, 0.00012598];
            AL_scaled=coeff(1)*(AL^10)+coeff(2)*(AL^9)+coeff(3)*(AL^8)+coeff(4)*(AL^7)+coeff(5)*(AL^6)+coeff(6)*(AL^5)+coeff(7)*(AL^4)+coeff(8)*(AL^3)+coeff(9)*(AL^2)+coeff(10)*(AL)+coeff(11);
        end        
    elseif label==7 %Lip_syn.
        if AL<0.005
            AL_scaled=AL;
        elseif AL>0.935
            AL_scaled=1;
        else
            coeff=[0, 0, 0, 0, -133.35, 408.21, -487.88, 287.27, -86.337, 12.329, 0.34625];
            AL_scaled=coeff(1)*(AL^10)+coeff(2)*(AL^9)+coeff(3)*(AL^8)+coeff(4)*(AL^7)+coeff(5)*(AL^6)+coeff(6)*(AL^5)+coeff(7)*(AL^4)+coeff(8)*(AL^3)+coeff(9)*(AL^2)+coeff(10)*(AL)+coeff(11);
        end
    end
end