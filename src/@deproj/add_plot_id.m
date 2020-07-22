function hts = add_plot_id( obj, ax )
%ADD_PLOT_ID Add the epicell ids to the specified plot axes.

    epicells = obj.epicells;
    n_objects = numel( epicells );
    hts = NaN( n_objects, 1 );
    
    for i = 1 :  n_objects
        o = epicells( i );
        center = double( o.center );
        hts(i ) = text( ax,  ...
            center(1), center(2), center(3) + 0.5, ...
            num2str( o.id ), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Clipping', 'on' );
    end
end

