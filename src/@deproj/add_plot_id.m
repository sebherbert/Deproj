function hts = add_plot_id( obj, ax )
%ADD_PLOT_ID Add the epicell ids to the specified plot axes.

    epicells = obj.epicells;
    n_objects = numel( epicells );
    hts = NaN( n_objects, 1 );
    
    for i = 1 :  n_objects
        o = epicells( i );
        hts(i ) = text( ax,  ...
            o.center(1), o.center(2), o.center(3) + 0.5, ...
            num2str( o.id ), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle' );
    end
end

