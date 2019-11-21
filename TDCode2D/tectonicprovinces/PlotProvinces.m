dir = 'C:\Research\tstar\TDCode2D\tectonicprovinces';

bnames={'Bound1'
        'Bound2'
        'Bound3'
        'Bound4'
        'Bound5'
        'Bound6'
        'Bound7'
        'Bound8'
        'Bound9'
        'Bound10'
        'Bound11'
        'Bound12'
        'Bound13'
        'Bound14'
        'Bound15'
        'Bound16'
        'Bound17'
        'Bound18'};
    
for knt=1:length(bnames)  
    
    B      = load([dir bnames{knt}]);
    [X, Y] = mfwdtran(mstruct,B(:,1),B(:,2));
    plot(X,Y,'k-','LineWidth',2);
    
end

