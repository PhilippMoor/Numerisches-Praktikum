function weight = Compute_weights_log()
    f=@(x) log(x);
    
    for k=1:45
        weight(k)=int(f,0,1);
    end
    
end