function [ ts, error, res ] = solve_ts_matrix( ts_vec, ind_1, ind_2, n_sta )
%solve_ts_matrix. weights is a vector of errors. can be used for dt
    
    A = zeros(length(ts_vec + 1), n_sta);
    
    for i = 1:length(ts_vec)
    
        A(i, ind_1(i)) = -1;
        A(i, ind_2(i)) = 1;
        
    end
    
    A(end + 1, :)   = 1;
    ts_vec(end + 1) = 0;%zero mean constraint
    
    ts = inv(A'*A)*A'*ts_vec';%solve
    
    for i = 1:(length(ts_vec)-1)
        
        res(i) = ts_vec(i) - (ts(ind_2(i)) - ts(ind_1(i)));
        
    end

    error = 0;%doesn't mean anything
    
end

