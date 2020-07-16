function hts = add_plot_distorsion_area( obj, ax )
%ADD_PLOT_DISTORSION_AREA Add the tissue plot colored with the error on
%cell area caused by the projection.
    
    epicells = obj.epicells;
    areas = vertcat( epicells.area );
    uncorr_areas = vertcat( epicells.uncorrected_area );
    err = 1 - uncorr_areas ./ areas;
    
    hts = add_plot_variable( obj, 100. * err, ax );
end

