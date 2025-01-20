#solve.jl

#formulate the problem at every solve
#pass in the initial condition, along with jacobians of the current reference trajectory 
#remove the bias term...
function update_prob(system, x_initial_k, all_Ad_k, all_Bd_k, P_k, unstable_directions_k, reference_traj_k, constraint_type)
    
    #state trajectory
    X = Convex.Variable(nx,N_h)
        
    #control trajectory 
    U = Convex.Variable(nu, N_h-1)
    
    #initial state constraint
    cons = Constraint[X[:,1] == x_initial_k]


    #dynamics constraint (first-order Taylor Expansion )
    for k=1:(N_h-1)
        #these are the dynamics for delta x
        #Linear affine term in terms of Δx and Δu (the integrator is there because the reference trajectory is not dynamically feasible)
        #(x̄+Δx) - (f(x̄,ū) + AΔx + BΔu) - bias
        push!(cons, zeros(6)== reference_traj_k[:,k+1] + X[:,k+1] - (dynamics_wcontrol_integrate(system, [reference_traj_k[:,k]; zeros(3)], Δt).u[end][1:6] + all_Ad_k[:,:,k]*X[:,k]+all_Bd_k[:,:,k]*U[1:3,k]))

    end

    if(constraint_type == "euc")

        #tube constraint
        #for k=2:(N_h-1) #this was there before and works
        #THIS WORKS and gives nice impluses. removes the noise. HIGH tube pose r
        for k=1:N_h 
           push!(cons, norm(X[1:3,k],2) <= tube_pose_r)
        
           push!(cons, norm(X[4:6,k],2) <= tube_vel_r)
            
        end

    else
    
        #cost to go constraint (Working, most current)
        for k=2:N_h
            
        costtogo = P_k[:,:,k]
            
        costtogo_hermitian = (costtogo + costtogo')/2
            
        #relaxing bc of the manifold constraint
        ctg_constraint = quadform(X[:,k], costtogo_hermitian) <= epsilon
            
        push!(cons, ctg_constraint)
                    
        end

    end

    #Contingency-constraint (half-space) (inequality)
    # for k = 2:N_h
                
    #     #manifold goes to the right
    #     manifold_constraint = X[:,k]'*unstable_directions_k[:,k]    

    #     #used to just be >, but there was a warning about it 
    #     push!(cons, manifold_constraint >= a)

    # end

    #strict equality
    for k = 2:N_h
                
        #manifold goes to the right
        manifold_constraint = X[:,k]'*unstable_directions_k[:,k]    

        #used to just be >, but there was a warning about it 
        push!(cons, manifold_constraint >= a)

    end

    
    return cons, X, U
        
end




#create the optimization problem
function solve_opt(cons, X, U, N)
        
    #the cost function (L1 norm on the cost and the tracking error)
    obj = 0

    for k=1:N-1
        
        obj += norm(U[:,k], 1)

    end            
    
    prob = minimize(obj, cons);
    
    #solve the problem with Mosek
    #silent_solver depracated, used silent instead 
    Convex.solve!(prob, ()-> Mosek.Optimizer(), silent= true);

    #obtain the state and control trajectories 
    Xm = X.value;
    Um = U.value;

    #display(prob.status)
    
    return Xm, Um
    
end