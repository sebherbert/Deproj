function hts = add_plot_curvature_gauss( obj, ax )
%ADD_PLOT_CURVATURE_GAUSS Add the tissue plot colored with the Gaussian curvature.
    
    epicells = obj.epicells;
    curvs = vertcat( epicells.curvatures );
    gauss_curv = curvs( :, 2 );
    
    hts = add_plot_variable( obj, gauss_curv, ax );
end

