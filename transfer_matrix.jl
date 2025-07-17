### Basic definitions of functions for 4d cosmological Regge calculus

function id(j0::Float64, j1::Float64, k::Float64)
    """
    Identity operator for insertion in the transfer matrix.
    """
    return 1.0
end

function height(j0::Float64, j1::Float64, k::Float64)
    """
    Height of 4-frustum as a function of the square areas j0 and j1 and the trapezoid area k.
    """
    res_squared = ((j0 - j1)^2 + 8*k^2)/(2 * (sqrt(j0) + sqrt(j1))^2)
    res = sqrt(complex(res_squared))
    return real(res)
end

function height_sq(j0::Float64, j1::Float64, k::Float64)
    """
    Squared height of 4-frustum as a function of the square areas j0 and j1 and the trapezoid area k.
    """
    res = ((j0 - j1)^2 + 8*k^2)/(2 * (sqrt(j0) + sqrt(j1))^2)
    return res
end

function four_volume(j0::Float64, j1::Float64, k::Float64)
    """
    4-Volume of 4-frustum as a function of the square areas j0 and j1 and the trapezoid area k.
    """
    res_squared = (j0 + j1)^2 / 32 * ((j0 - j1)^2 + 8*k^2)
    res = sqrt(complex(res_squared))
    return real(res)
end

function deficit_angle_trap(j0::Float64, j1::Float64, k::Float64)
    """
    Deficit angle of 4-frustum located at trapezoids as a function of the square areas j0 and j1 and the trapezoid area k.
    """
    return (pi/2 - acos((j0 - j1)^2 / (16*k^2 + (j0 - j1)^2)))
end

function deficit_angle_square(j0::Float64, j1::Float64, k::Float64)
    """
    Deficit angle of 4-frustum located at squares as a function of the square areas j0 and j1 and the trapezoid area k.
    """
    return real(asinh((j1 - j0) / (sqrt(complex(16*k^2 + (j1 - j0)^2)))))
end

function S_III(j0::Float64, j1::Float64, k::Float64, Λ::Float64) 
    """
    Single 4-frustum Regge action in the causally regular sector (Sector III) as a function of the square areas j0 and j1 and the trapezoid area k.
    Λ: cosmological constant
    """
    res = 6 * (j0 - j1) * asinh((j1 - j0) / (sqrt(complex(16*k^2 + (j1 - j0)^2)))) + 12 * k * (pi/2 - acos((j0 - j1)^2 / (16*k^2 + (j0 - j1)^2))) - Λ * four_volume(j0, j1, k)
    return real(res)
end

### Solving the equations of motion

function dS_IIIdk(j0::Float64, j1::Float64, k::Float64, Λ::Float64)
    """
    Derivative of 4-frustum Regge action S_III with respect to trapezoid area for equations of motion.
    Λ: cosmological constant
    """
    res = 12 * (pi/2 - acos((j0 - j1)^2 / (16*k^2 + (j0 - j1)^2))) - Λ * ((j0 + j1)^2 * k) / (4 * four_volume(j0, j1, k))
    return res
end

function dS_IIIdj1(j0::Float64, j1::Float64, j2::Float64, k0::Float64, k1::Float64, Λ::Float64)
    """
    Derivative of 4-frustum Regge action S_III with respect to square area j1 for equations of motion. This involves two 4-frusta.
    Λ: cosmological constant
    """
    dS0dj1 = -6  * asinh((j1 - j0) / (sqrt(complex(16 * k0^2 + (j1 - j0)^2)))) - Λ * (j1 * (j1 - j0) + 4 * k0^2) / (2 * sqrt(2) * sqrt((j0 - j1)^2 + 8 * k0^2))
    dS1dj1 = -6  * asinh((j1 - j2) / (sqrt(complex(16 * k1^2 + (j1 - j2)^2)))) - Λ * (j1 * (j1 - j2) + 4 * k1^2) / (2 * sqrt(2) * sqrt((j1 - j2)^2 + 8 * k1^2))
    return real(dS0dj1 + dS1dj1)
end

function eom_nslice(j0::Float64, jf::Float64, Λ::Float64, γ::Float64, V::Int64)
    """
    Generate the area Regge equations of motion.
    Input:
        j0: initial square area.
        jf: final square area.
        Λ: cosmological constant.
        γ: Barbero-Immirzi parameter.
        V: number of building blocks in temporal direction (V) and spatial direction (V^3).
    Output:
        eom_vector: function that takes a vector of spins (k0, j1, k1, j2, ..., k_{V-1}) and maps it to equations of motion. 
    Formatted for input of NLsolve.jl.
    """
    function eom_vector(blk_spins::Vector{Float64})
        """
        Return array of equations of motion for input spins (k0, j1, k1, j2, ..., k_{V-1}).
        """
        eom_vec = []
        for n in 1:Int(2 * V - 1)
            if n == 1 && V > 1
                push!(eom_vec, dS_IIIdk(γ * j0, γ * blk_spins[2], blk_spins[1], Λ))
            elseif n == 1 && V == 1
                push!(eom_vec, dS_IIIdk(γ * j0, γ * jf, blk_spins[1], Λ))
            elseif n == 2 && V > 2
                push!(eom_vec, dS_IIIdj1(γ * j0, γ * blk_spins[2], γ * blk_spins[4], blk_spins[1], blk_spins[3], Λ))
            elseif n == 2 && V == 2
                push!(eom_vec, dS_IIIdj1(γ * j0, γ * blk_spins[2], jf, blk_spins[1], blk_spins[3], Λ))
            elseif n ==  (2 * V - 2)
                push!(eom_vec, dS_IIIdj1(γ * blk_spins[2 * (V - 2)], γ * blk_spins[2 * (V - 1)], γ * jf, blk_spins[2 * V - 3], blk_spins[2 * V - 1], Λ))
            elseif n == (2 * V - 1)
                push!(eom_vec, dS_IIIdk(γ * blk_spins[2 * (V - 1)], γ * jf, blk_spins[2 * V - 1], Λ))
            elseif isodd(n)
                push!(eom_vec, dS_IIIdk(γ * blk_spins[n - 1], γ * blk_spins[n + 1], blk_spins[n], Λ))
            elseif iseven(n)
                push!(eom_vec, dS_IIIdj1(γ * blk_spins[n - 2], γ * blk_spins[n], γ * blk_spins[n + 2], blk_spins[n - 1], blk_spins[n + 1], Λ))
            end
        end
        return eom_vec
    end
    return eom_vector
end

function generate_initguess(j0::Float64, jf::Float64, k_guess::Float64, Λ::Float64, γ::Float64, V::Int64)
    """
    Generate inital guess for the solutions to the equations of motion.
    Input:
        j0: initial square area.
        jf: final square area.
        k_guess: initial guess for the trapezoid area solution for a single 4-frustum.
        Λ: cosmological constant.
        γ: Barbero-Immirzi parameter.
        V: number of building blocks in temporal direction (V) and spatial direction (V^3).
    """
    k_sol = nlsolve(eom_1slice(j0, jf, Λ, γ), [k_guess], iterations=10000, xtol=1e-10, ftol=1e-10)
    if converged(k_sol) == true
        k = k_sol.zero[1]
    else
        return "Obtaining the 1-slice k solution failed. Try different k_guess!"
    end
    dj = (jf - j0) / V
    init_vec = Float64[]
    for n in 1:(2 * V - 1)
        if isodd(n)
            push!(init_vec, k/V)
        elseif iseven(n)
            #push!(init_vec, j_evol_class(j0, jf, Int(n / 2), V)) using here the continuum equation to guess
            push!(init_vec, j0 + (n / 2) * dj)
        end
    end
    return init_vec
end

function solve_eom(j0::Float64, jf::Float64, Λ::Float64, γ::Float64, V::Int64, k_guess::Float64, iterations::Int64=10000, xtol::Float64=1e-10, ftol::Float64=1e-10)
    """
    Solve the equations of the area Regge equations of motion via nsolve method of NLsolve.
    Input:
        j0: initial square area.
        jf: final square area.
        Λ: cosmological constant.
        γ: Barbero-Immirzi parameter.
        V: number of building blocks in temporal direction (V) and spatial direction (V^3).
        k_guess: initial guess for the trapezoid area solution for a single 4-frustum.
        iterations: maximum number of iterations for the solver.
        xtol: tolerance for the solution vector.
        ftol: tolerance for the function value.
    Output:
        js: solution of the square areas.
        ks: solution of the trapezoid areas.
        Hs: heights of the 4-frusta.
        prop_time: accumulated height of 4-frusta.
    """
    discr_sols = nlsolve(eom_nslice(j0, jf, Λ, γ, V), generate_initguess(j0, jf, k_guess, Λ, γ, V), iterations=iterations, xtol=xtol, ftol=ftol).zero
    js = vcat([j0], discr_sols[2:2:end], [jf])
    ks = discr_sols[1:2:end]
    Hs = [height(js[n], js[n+1], ks[n]) for n in 1:V]
    prop_time = vcat([0.0], accumulate(+, Hs))
    return (js, ks, Hs, prop_time)
end

### The effective spin-foam path integral

function μ_toy(j0::Float64, j1::Float64, k::Float64, α::Float64)
    """
    Local measure for the path integral as a function of (j0, j1, k).
    α: exponent of the denominator, controlling fall-off behavior in k and j. For absolute convergence of path integral with k insertion \alpha > 3/2.
    """
    return k  / ((j0 - j1)^2 + 8 * k^2)^(α)
end

function periodicity(j0, j1, Λ)
    """
    Periodicity of oscillations of Regge action from linearization.
    Use to determine necessary refinement of the spectrum of k to resolve oscillations and stationary points.
    """
    return 4*pi / (Λ * (j0 +  j1))
end

function linear_k_sol(j0, j1, Λ)
    """
    Solution of the linearized Regge equations for the trapezoid area k. 
    Use to estimate kmax for trans_mat function.
    """
    return abs(j0 - j1) * real(sqrt(complex((3 + Λ * (j0 + j1) / 8) / (2 * Λ * (j0 + j1))))) / 2
end

function trans_mat(O::Function, N::Int64, γ::Float64, Λ::Float64, rtol::Float64, V::Int64=1, ref_sl::Float64=1.0, ref_tl::Float64=1.0, α::Float64=2.5)
    """
    Compute the transfer matrix, with indices corresponding to boundary areas. For every component, perform a summation over the timelike trapezoid areas via Wynn's epsilon algorithm.
    Input:
        O: operator as a function of (j0, j1, k) to be inserted into the path integral.
        jmax: determines the range of square areas j0 and j1, which are in the range [0, jmax].
        γ: Barbero-Immirzi parameter.
        Λ: cosmological constant.
        rtol: relative tolerance for convergence of Wynn's epsilon algorithm.
        V: number of building blocks in temporal direction (V) and spatial direction (V^3).
        ref_sl: refinement factor for spectrum of spacelike area.
        ref_tl: refinement factor for spectrum of timelike area.
        α: parameter of the measure.
    Output:
        mat: transfer matrix with general operator insertion O.
    N.B.: 
        - A refinement of the spectrum by 1/V^2 is already included.
        - Wynn's algorithm is cut off at kmax = 
        - The measure is defined spatially 'non-local' and takes the full boundary spinss j0 and j1 into account.
    """
    mat = zeros(ComplexF64, N, N)
    k_spec = periodicity()
    for i0 in 1:N
        for i1 in 1:j0
            kmax = (ref_tl / V^2 * 10 * abs(j0 - j1) + 50)
            kvals = (ref_tl / V^2 * 0.5):(ref_tl / V^2 * 0.5):kmax
            amplitudes = [O(ref_sl * γ * 0.5 * j0 / V^2, ref_sl * γ * 0.5 * j1 / V^2, k) * μ_toy(ref_sl * γ * 0.5 * j0, ref_sl * γ * 0.5 * j1, k, α) *  exp(V^3 * im * S_III(ref_sl * γ * 0.5 * j0 / V^2, ref_sl * γ * 0.5 * j1 / V^2, k, Λ)) for k in kvals]
            mat_element = wynn(accumulate(+, amplitudes)[end - 50:end], rtol)[1]
            mat[j0, j1] = mat_element
            mat[j1, j0] = mat_element
        end
    end
    return mat
end