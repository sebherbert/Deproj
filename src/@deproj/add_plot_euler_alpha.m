function hts = add_plot_euler_alpha( obj, ax )
%ADD_PLOT_EULER_ALPHA Add the tissue plot colored with the 1st euler angle.
    
    epicells = obj.epicells;
    angles = vertcat( epicells.euler_angles );
    alphas = rad2deg( angles( :, 1 ) );
    
    % Wrap back to 0 - 180º.
    neg_alphas = alphas < 0;
    alphas( neg_alphas ) = 180 + alphas( neg_alphas );

    hts = add_plot_variable( obj, alphas, ax );
end

