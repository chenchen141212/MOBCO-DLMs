
function eta_=geteta(eta_,S)
    if rand<=0.5
        eta=eta_+0.1*tan(pi*(rand-0.5)); % cauchy distribution
    else
        eta=normrnd(eta_,0.1);
    end    
    eta_=(eta_*S+eta)/(S+1);

end

    