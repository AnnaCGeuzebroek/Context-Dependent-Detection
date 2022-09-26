% check for button down

[xblah,yblah,buttons] = GetMouse(whichScreen); 
if any(buttons)
    if ButtonDown==0
        numResp=numResp+1; 
        RespT(numResp)=GetSecs;
       % RespTcopy(numResp)=GetSecs;
        if buttons(1)
            RespLR(numResp)=1;
            disp('Button 1')
        elseif buttons(end)
            RespLR(numResp)=2;
            disp('Button 2')
        end
        if par.recordEEG, sendtrigger64(par.CD_BUTTONS(RespLR(numResp)),port,ioObj); end
    end    
    ButtonDown = 1;
else
    ButtonDown = 0;
end