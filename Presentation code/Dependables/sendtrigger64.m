function sendtrigger64(trg,port,ioObj)

io32(ioObj,port,trg);
WaitSecs(0.003);
io32(ioObj,port,0);

 