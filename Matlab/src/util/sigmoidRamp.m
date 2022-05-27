function k = sigmoidRamp(tzero,tpause,haltingiter)
    if haltingiter - tpause < tzero
        k = 0;
    else
        t = haltingiter - tzero - tpause;
        
        c = (haltingiter - tzero)/2;
        
        k = 1/(1+exp(-(10/c)*t + 10));
    end
end