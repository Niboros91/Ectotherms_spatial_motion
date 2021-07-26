function A_tot = compute_A_multi_16(N_trees,Pest_stages,stages, Adj,u,trap,trap_mortality)
    
    %Total number of states
    n = N_trees*(stages);
    
    indices = find(u>0); %Nodes that have a trap
   
    %Compute the A matrix for each area
    for i = 1:N_trees
        
        A_sing{i} = [-Pest_stages(1,i).growth-Pest_stages(1,i).death 0 0 0 0 0 0  Pest_stages(8,i).birth 0 0 0;... %Egg (1)
            Pest_stages(1,i).growth -Pest_stages(2,i).growth-Pest_stages(2,i).death  0 0 0 0 0 0 0 0 0;... %L1 (2)
            0 Pest_stages(2,i).growth -Pest_stages(3,i).growth-Pest_stages(3,i).death  0 0 0 0 0 0 0 0;... %L2 (3)
            0 0 Pest_stages(3,i).growth -Pest_stages(4,i).growth-Pest_stages(4,i).death  0 0 0 0 0 0 0;... %L3 (4)
            0 0 0  Pest_stages(4,i).growth -Pest_stages(5,i).growth-Pest_stages(5,i).death 0 0 0 0 0 0;...%P (5)
            0 0 0 0  Pest_stages(5,i).sex_ratio*Pest_stages(6,i).growth -Pest_stages(6,i).growth-Pest_stages(6,i).death 0 0 0 0 0;...%AM (6)
            0 0 0 0  (1-Pest_stages(6,i).sex_ratio)*Pest_stages(7,i).growth 0 -(Pest_stages(7,i).mate-Pest_stages(7,i).death) 0 0 0 0;...%NMF (7)
            0 0 0 0 0 0 Pest_stages(7,i).mate-Pest_stages(7,i).death -Pest_stages(8,i).growth-Pest_stages(8,i).death 0 0 0;... %MF (8)
            0 0 0 0 0 0 0 0 0 0 0;... %Traps AM (9)
            0 0 0 0 0 0 0 0 0 0 0;... %Traps NMF (10)
            0 0 0 0 0 0 0 0 0 0 0]; %Traps MF (11)


        %We have to use u and adj

        
    end
    
    for i = indices'
       if ~isempty(indices)
        indices_2 = find(Adj(i,:)>0);
            if ~isempty(indices_2)
                for j = indices_2
                    % We add the loses to the adult stages
                    A_sing{j}(6,6) = A_sing{j}(6,6) - trap;
                    A_sing{j}(7,7) = A_sing{j}(7,7) - trap ;
                    A_sing{j}(8,8) = A_sing{j}(8,8) - trap;
                end
            end
       end
       % Mortality due to the trap
       A_sing{i}(6,6) = A_sing{i}(6,6) - trap_mortality;
       A_sing{i}(7,7) = A_sing{i}(7,7) - trap_mortality;
       A_sing{i}(8,8) = A_sing{i}(8,8) - trap_mortality; 
       
       % Increase in the catching of the traps
       A_sing{i}(9,6) =  trap_mortality;
       A_sing{i}(10,7) =  trap_mortality;
       A_sing{i}(11,8) =  trap_mortality; 
    end

    %Create a block diagonal matrix based on each areas A matrix
    A = blkdiag(A_sing{:});
    
    %% Adding the non diagonal elements
    
    %Dummy matrix
    
    A_jump = [0 0 0 0 0 0 0 0 0 0 0;...%Egg (1)
        0 0 0 0 0 0 0 0 0 0 0;... %L1 (2)
        0 0 0 0 0 0 0 0 0 0 0;... %L2 (3)
        0 0 0 0 0 0 0 0 0 0 0;... %L3 (4)
        0 0 0 0 0 0 0 0 0 0 0;... %P (5)
        0 0 0 0 0 1 0 0 0 0 0;... %AM (6)
        0 0 0 0 0 0 1 0 0 0 0;... %NMF (7)
        0 0 0 0 0 0 0 1 0 0 0;... %MF (8)
        0 0 0 0 0 0 0 0 0 0 0;... %Traps AM (9)
        0 0 0 0 0 0 0 0 0 0 0;... %Traps NMF (10)
        0 0 0 0 0 0 0 0 0 0 0]; %Traps MF (11) 
    
    A_jump_tot = zeros(n,n);
    
    %Based on the wind and the distribution we add the non diagonal
    %elements
     for i = indices'
        indices_2 = find(Adj(i,:)>0);
        for k = indices_2
            
            A_jump_temp = A_jump.*trap;
            temp_Adj = zeros(N_trees,N_trees);
            temp_Adj(i,k) = 1 ; % We insert it at the jump interaction
            A_jump_ind = kron(temp_Adj,A_jump_temp); %We expand to the overall system
            
            %We add to the big matrix
            A_jump_tot = A_jump_tot + A_jump_ind;
            
        end
     end

    %The final matrix is the sum of the block diagonal matrix and the interactions matrix 
    A_tot = A + A_jump_tot; 
  
    
end