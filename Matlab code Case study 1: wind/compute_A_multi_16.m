function A_tot = compute_A_multi_16(N_trees,Pest_stages,stages, Adj_w,w)
    
    %Total number of states
    n = N_trees*stages;
 
    %Check the directions of the wind
    sign_w=sign(w);
    
    %Compute the A matrix for each area
    for i = 1:N_trees
        
        A_sing{i} = [-Pest_stages(1,i).growth-Pest_stages(1,i).death 0 0 0 0 0 0  Pest_stages(8,i).birth;... %Egg (1)
            Pest_stages(1,i).growth -Pest_stages(2,i).growth-Pest_stages(2,i).death  0 0 0 0 0 0;... %L1 (2)
            0 Pest_stages(2,i).growth -Pest_stages(3,i).growth-Pest_stages(3,i).death  0 0 0 0 0;... %L2 (3)
            0 0 Pest_stages(3,i).growth -Pest_stages(4,i).growth-Pest_stages(4,i).death  0 0 0 0;... %L3 (4)
            0 0 0  Pest_stages(4,i).growth -Pest_stages(5,i).growth-Pest_stages(5,i).death 0 0 0;...%P (5)
            0 0 0 0  Pest_stages(5,i).sex_ratio*Pest_stages(6,i).growth -Pest_stages(6,i).growth-Pest_stages(6,i).death 0 0;...%AM (6)
            0 0 0 0  (1-Pest_stages(6,i).sex_ratio)*Pest_stages(7,i).growth 0 -(Pest_stages(7,i).mate-Pest_stages(7,i).death) 0;...%NMF (7)
            0 0 0 0 0 0 Pest_stages(7,i).mate-Pest_stages(7,i).death -Pest_stages(8,i).growth-Pest_stages(8,i).death]; %MF (8)

        %Based on the wind (and the position of the area) the area loses
        %insects or not
            temp_jump =0;
            for k =1:N_trees
            
                if  Adj_w(i,k)*sign_w(1) == 2
                    
                        temp_jump = temp_jump  - abs(w(1));
                end
                if  Adj_w(i,k)*sign_w(2) == 1
                        temp_jump = temp_jump  - abs(w(2));
                end
                
            end
            
            % We add the loses to the adult stages
            A_sing{i}(6,6) = A_sing{i}(6,6) + temp_jump;
            A_sing{i}(7,7) = A_sing{i}(7,7) + temp_jump ;
            A_sing{i}(8,8) = A_sing{i}(8,8) + temp_jump;
        
    end

    %Create a block diagonal matrix based on each areas A matrix
    A = blkdiag(A_sing{:});
    
    %% Adding the non diagonal elements
    
    %Dummy matrix
    
    A_jump = [0 0 0 0 0 0 0 0;...%Egg (1)
        0 0 0 0 0 0 0 0;... %L1 (2)
        0 0 0 0 0 0 0 0;... %L2 (3)
        0 0 0 0 0 0 0 0;... %L3 (4)
        0 0 0 0 0 0 0 0;... %P (5)
        0 0 0 0 0 1 0 0;... %AM (6)
        0 0 0 0 0 0 1 0;... %NMF (7)
        0 0 0 0 0 0 0 1]; %MF (8)
    
    A_jump_tot = zeros(n,n);
    
    %Based on the wind and the distribution we add the non diagonal
    %elements
     for i = 1:N_trees 
        for k = 1:N_trees
            temp_jump = 0;
            
            if  Adj_w(i,k)*sign_w(1) == 2
                    
                     temp_jump =  abs(w(1));
            end
            if  Adj_w(i,k)*sign_w(2) == 1
                     temp_jump =  abs(w(2));
            end
            
            A_jump_temp = A_jump.*temp_jump;
            
            
            temp_Adj = zeros(N_trees,N_trees);
            temp_Adj(k,i) = 1 ; % We insert it at the jump interaction
            A_jump_ind = kron(temp_Adj,A_jump_temp); %We expand to the overall system
            
            %We add to the big matrix
            A_jump_tot = A_jump_tot + A_jump_ind;
            
        end
     end

    %The final matrix is the sum of the block diagonal matrix and the interactions matrix 
    A_tot = A + A_jump_tot; 
  
    
end