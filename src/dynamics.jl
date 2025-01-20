#three-body problem dynamics
#these functions are dependent on a three body system that is defined by parameters.jl

function effective_potential(system, X)
    x = X[1]
    y = X[2]
    z = X[3]

    r1 = sqrt((x+system.μ2)^2 + y^2 + z^2) 
    r2 = sqrt((x-system.μ1)^2 + y^2 + z^2)
    #assuming m3 is unit mass
   
    U = (-system.μ1/r1)-(system.μ2/r2)-0.5*(x^2+y^2)
   
    return U
    
end


#three body problem dynamics with no control 
#these dynamics are in the cr3bp scaled units 
function three_body_prob_dynamics(system, x)
    
    q = zeros(eltype(x),3)
    v = zeros(eltype(x),3)
    q = x[1:3]
    v = x[4:6]

    ẋ = zeros(eltype(x),6)
    U_q = zeros(eltype(x),3)
    
    U_q = (ForwardDiff.gradient(_x -> effective_potential(system, _x), q))
        
    ẋ[1:3] = v
    ẋ[4] = 2*v[2] - U_q[1] 
    ẋ[5] = -2*v[1] - U_q[2]
    ẋ[6] = -U_q[3]

    return ẋ
end


#three-body dynamics with control 
#x[1:6] -> state
#x[7:9] -> control inputs
#these dynamics are in the cr3bp scaled units 
function three_body_prob_dynamics_wcontrol(system, x)
    
    q = zeros(eltype(x),3)
    v = zeros(eltype(x),3)
    q = x[1:3]
    v = x[4:6]
    u = x[7:9]

    ẋ = zeros(eltype(x),9)
    
    U_q = zeros(eltype(x),3)
    
    U_q = (ForwardDiff.gradient(_x -> effective_potential(system, _x), q))
        
    ẋ[1:3] = v
    ẋ[4] = 2*v[2] - U_q[1] + u[1] 
    ẋ[5] = -2*v[1] - U_q[2] + u[2]
    ẋ[6] = -U_q[3] + u[3]
    
    #zero order hold on the controls 
    ẋ[7:9] = zeros(3)

    return ẋ
end


#these dynamics the position is in km and time in days 
function three_body_prob_dynamics_scaled(system, x)

    q_original = zeros(eltype(x),3)
    v_original = zeros(eltype(x),3)
    
    q_original = x[1:3]/system.position_scale 
    v_original = x[4:6]/system.velocity_scale
    
    x_original = [q_original; v_original]

    #original is in the CR3BP units
    ẋ_original = zeros(eltype(x),6)

    #calculate the original xdot (no scaling)

    ẋ_original = three_body_prob_dynamics(system, x_original)

    #then scale the output
    v_scaled = ẋ_original[1:3]*system.velocity_scale
    
    a_scaled = ẋ_original[4:6]*system.acceleration_scale

    ẋ_scaled = [v_scaled; a_scaled]

    return ẋ_scaled

end


#inputs are in the custom scaled units
function three_body_prob_dynamics_wcontrol_scaled(system, x)

    q_original = zeros(eltype(x),3)
    v_original = zeros(eltype(x),3)
    u_original = zeros(eltype(x),3)
    
    q_original = x[1:3]/system.position_scale 
    v_original = x[4:6]/system.velocity_scale
    u_original = x[7:9]/system.acceleration_scale
    
    x_original = [q_original; v_original; u_original]

    ẋ_original = zeros(eltype(x),9)

    #calculate the original xdot (no scaling)

    #xoriginal is in the CR3BP units
    ẋ_original = three_body_prob_dynamics_wcontrol(system, x_original)

    #then scale the output
    v_scaled = ẋ_original[1:3]*system.velocity_scale
    
    a_scaled = ẋ_original[4:6]*system.acceleration_scale

    u_scaled = ẋ_original[7:9]*system.acceleration_scale

    ẋ_scaled = [v_scaled; a_scaled; u_scaled]

    return ẋ_scaled

end


#inputs are in the custom scaled units

function three_body_prob_dynamics_wcontrol_scaled_RK4(system, x, u)

    q_original = zeros(eltype(x),3)
    v_original = zeros(eltype(x),3)
    u_original = zeros(eltype(x),3)
    
    q_original = x[1:3]/system.position_scale 
    v_original = x[4:6]/system.velocity_scale
    u_original = u/system.acceleration_scale
    
    x_original = [q_original; v_original; u_original]

    ẋ_original = zeros(eltype(x),9)

    #calculate the original xdot (no scaling)

    #xoriginal is in the CR3BP units
    ẋ_original = three_body_prob_dynamics_wcontrol(system, x_original)

    #then scale the output
    v_scaled = ẋ_original[1:3]*system.velocity_scale
    
    a_scaled = ẋ_original[4:6]*system.acceleration_scale

    u_scaled = ẋ_original[7:9]*system.acceleration_scale

    ẋ_scaled = [v_scaled; a_scaled]

    return ẋ_scaled

end