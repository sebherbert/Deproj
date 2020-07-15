function hts = add_plot_curvature_mean( obj, ax )
%ADD_PLOT_CURVATURE_MEAN Add the tissue plot colored with the mean curvature.
    
    epicells = obj.epicells;
    curvs = vertcat( epicells.curvatures );
    mean_curv = curvs( :, 1 );
    
    hts = add_plot_variable( obj, mean_curv, ax );
end

