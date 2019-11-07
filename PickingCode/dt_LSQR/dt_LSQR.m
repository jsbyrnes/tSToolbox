function [ dt_solved, res ] = dt_LSQR( Traces )

    n_traces = length(Traces);
    n_pairs = n_traces^2;
    
    for i = 1:(n_pairs)

        ind1(i) = mod(i-1, n_traces) + 1;
        ind2(i) = floor((i-1)/n_traces) + 1;
        dt(i)   = pair_measurement( Traces(ind1(i)), Traces(ind2(i)));
                
    end

    [dt_solved, ~, res]   = solve_ts_matrix(dt, ind1, ind2, n_traces);
    
end

