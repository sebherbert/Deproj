function hts = add_plot_distorsion_perimeter( obj, ax )
%ADD_PLOT_DISTORSION_PERIMETER Add the tissue plot colored with the error on
%cell perimeter caused by the projection.
    
    epicells = obj.epicells;
    perims = vertcat( epicells.perimeter );
    uncorr_perims = vertcat( epicells.uncorrected_perimeter );
    err = 1 - uncorr_perims ./ perims;
    
    hts = add_plot_variable( obj, 100. * err, ax );
end

