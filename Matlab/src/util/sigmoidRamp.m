function k = sigmoidRamp(tzero,tpause,haltingiter)
    if haltingiter - tpause < tzero
        k = 0;
    else
        t = haltingiter - tzero - tpause;
        
        c = (haltingiter - tzero)/2;
        
        k = 1/(1+exp(-(6/c)*t + 6));
    end
end