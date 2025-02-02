#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# This julia file provides an explicit implementation of Wynn's epsilon algorithm to accelerate sequence convergence. Includes method for more than one summation variable #
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

function acc_triang(mat::Matrix{ComplexF64}, N::Int64)

    res = zeros(ComplexF64, N)

    for n in 1:N

        part_sum = 0.0 + 0.0im
        
        for i in 1:n

            for j in 1:(n-i+1)
                
                part_sum += mat[i,j]

            end

        end

        res[n] = part_sum

    end

    res

end

function acc_rect(mat::Matrix{ComplexF64}, N::Int64)

    res = zeros(ComplexF64, N)

    for n in 1:N

        part_sum = 0.0 + 0.0im
        
        for i in 1:n

            for j in 1:n
                
                part_sum += mat[i,j]

            end

        end

        res[n] = part_sum

    end

    res

end

function wynn(seq::Vector{ComplexF64}) 

    n = length(seq)
    e = zeros(ComplexF64, n, n+1)

    for i in 1:n

        e[i, 2] = seq[i]

    end
    
    for k in 3:(n+1)
        for i in 1:(n+2-k)
            
            e[i, k] = e[i+1, k-2] + 1 / (e[i+1, k-1] - e[i, k-1])
        
        end
    end
    
    if iseven(n)

        res = e[1,n]

    else
        
        res = e[1,n-1]
        
    end
        
    res
    
end