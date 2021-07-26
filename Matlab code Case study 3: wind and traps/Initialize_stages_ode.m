function Pest_stages = Initialize_stages_ode(Birth,Death,Growth,S_R,R_mate,R_remate,Pest_stages,N_trees,N_stages)
    %Input =  Different rates at each instant
    %Output =  Parameters associated to each stage
    for j=1:N_trees
        for i=1:N_stages %At this point all stages have the same parameters
            Pest_stages(i,j).birth = Birth(j);
            Pest_stages(i,j).death = Death(j);
            Pest_stages(i,j).growth = Growth(j);
            Pest_stages(i,j).sex_ratio = S_R(j);
            Pest_stages(i,j).mate = R_mate(j);
            Pest_stages(i,j).remate = R_remate(j);

        end
    end

end