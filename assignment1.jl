using LinearAlgebra

function choleski_decomposition(A::Array{Float64},b::Array{Float64})
    n = size(A,1)
    L = zeros(n,n)
    y = zeros(n,1)
    x = zeros(n,1)
    #decompose L

    for j = 1:n

        sumⱼ = 0
        for q = 1:j-1
            sumⱼ = sumⱼ + L[j,q]^2
        end
        L[j,j] = sqrt(A[j,j] - sumⱼ)

        sumᵢ = 0
        for i = (j+1):n

            for k = 1:j-1
                if k == 1
                    sumᵢ = L[i,k]*L[j,k]
                else
                    sumᵢ = sumᵢ + L[i,k]*L[j,k]
                end
            end
            L[i,j] = (A[i,j] - sumᵢ) / L[j,j]
        end
    end

    #get y
    for i = 1:n
        sumᵧ = 0
        for j = 1:i-1
            sumᵧ = sumᵧ + L[i,j]*y[j]
        end
        y[i] = (b[i] - sumᵧ) / L[i,i]
    end
    #calculate x
    for i = n:-1:1
        sumₓ = 0
        for j = (i+1):n
            sumₓ = sumₓ + L[j,i]*x[j]
        end
        x[i] = (y[i] - sumₓ) / L[i,i]
    end

    return L,y,x
end
