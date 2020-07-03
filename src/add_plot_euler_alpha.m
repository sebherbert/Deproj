function hts = add_plot_euler_alpha( ax, epicells )
%ADD_PLOT_EULER_ALPHA Add the tissue plot colored with the 1st euler angle.
    
    angles = vertcat( epicells.euler_angles );
    alphas = rad2deg( angles( :, 1 ) );
    
    % Wrap back to 0 - 180º.
    neg_alphas = alphas < 0;
    alphas( neg_alphas ) = 180 + alphas( neg_alphas );

    hts = add_plot_variable( ax, { epicells.boundary }, alphas );
end

