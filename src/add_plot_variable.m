function hts = add_plot_variable( ax, boundaries, values )
%ADD_PLOT_VARIABLE Plots the boundaries as patches, colored by the specified values.
    
    n_objects = numel( boundaries );
    hts = NaN( n_objects, 1 );
    
    for i = 1 :  n_objects

        p = boundaries{ i };
        val = values( i );
        hts(i ) = patch( ax, p(:,1), p(:,2), p(:,3), val, ...
            'LineWidth', 2 );
    end
end
