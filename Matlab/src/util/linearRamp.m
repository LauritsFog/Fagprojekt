function k = linearRamp(tzero,tpause,haltingiter)
    if haltingiter - tpause < tzero
        k = 0;
    else
        k = (haltingiter - tzero - tpause)/(haltingiter - tzero);
    end
end