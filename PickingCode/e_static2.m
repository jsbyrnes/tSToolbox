function [ e_static ] = e_static( all_tS, dataE )


    Gevt = zeros(length(all_tS), length(unique(dataE)));
    %build row-by-row
    
    for k=1:length(all_tS)
        
       Gevt(k,:)=0;
       
       Gevt(k,dataE(k))=1;
    
    end

    e_static = G\all_tS;

end

