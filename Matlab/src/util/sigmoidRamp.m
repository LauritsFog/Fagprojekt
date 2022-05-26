function k = sigmoidRamp(tzero,tpause,haltingiter)
    if haltingiter - tpause < tzero
        k = 0;
    else
        t = haltingiter - tzero - tpause;
        
        k = 1/(1+exp(-(t - (haltingiter - tzero)/2)));
    end
end