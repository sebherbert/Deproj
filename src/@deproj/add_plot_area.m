function hts = add_plot_area( obj, ax )
%ADD_PLOT_AREA Add the tissue plot colored with the cell area.
    
    epicells = obj.epicells;
    areas = vertcat( epicells.area );
    
    hts = add_plot_variable( obj, areas, ax );
end

