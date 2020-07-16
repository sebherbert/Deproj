function hts = add_plot_perimeter( obj, ax )
%ADD_PLOT_PERIMETER Add the tissue plot colored with the cell perimeter.
    
    epicells = obj.epicells;
    perims = vertcat( epicells.perimeter );
    
    hts = add_plot_variable( obj, perims, ax );
end

