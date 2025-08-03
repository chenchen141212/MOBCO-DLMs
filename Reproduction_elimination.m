
function pop=Reproduction_elimination(rep,pop,nPoP,Sr)
   
    temppop = pop;

    % Sort temppop based on ReporPopInd in ascending order
    [~, sortIdx] = sort([temppop.ReporPopInd]);
    temppop = temppop(sortIdx);

    % Combine rep and temppop
    learnedpop = [rep; temppop];

    for i=1:nPoP
        pop(i).Position=learnedpop(i).Position;
        pop(i).StepSize=learnedpop(i).StepSize;
    end

    for i=1:Sr
        pop(nPoP+1-i).Position=pop(i).Position;
        pop(nPoP+1-i).StepSize=pop(i).StepSize;
    end   

end