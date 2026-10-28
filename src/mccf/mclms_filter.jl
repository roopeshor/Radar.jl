"""
    mclms_filter(X_r, X_r_prime, mclms_len)

# Arguments
- `X_r::AbstractVector{<:Number}`
- `X_r_prime::AbstractVector{<:Number}`
- `mclms_len::Integer`

# Returns
- `Y_r::AbstractVector`
"""
function mclms_filter(X_r::AbstractVector{<:Number}, X_r_prime::AbstractVector{<:Number}, mclms_len::Integer)
    N2 = length(X_r)
    # Normalized step size u
    u = 1.0 / sum(abs.(X_r) .^ 2)
    # Initial filter weights W1 = 0
    W = zeros(eltype(X_r), mclms_len)

    # Pass 1: Adaptive coefficient training
    for j in 1:N2
        if j < mclms_len
            X_vec = vcat(zeros(eltype(X_r), mclms_len - j), X_r[1:j])
        else
            X_vec = X_r[j - mclms_len + 1 : j]
        end
        Y = dot(W, X_vec)
        e_j = X_r_prime[j] - Y
        W .= 0.98 .* W .+ u * e_j * conj.(X_vec)          # Weight update
    end

    # Pass 2: Re-filtering with updated weights W
    # Y_r = zeros(1, N2);
    # for j = 1:N2
    #   if j < mclms_len
    #       X_vec = [zeros(1, mclms_len - j), X_r(1:j)];
    #   else
    #       X_vec = X_r(j - mclms_len + 1 : j);
    #   end
    #   Y_r(j) = W * X_vec.';
    # end

    # 0.6ms/60% faster than above
    b = reverse(W)
    Y_r = conv(b, X_r)[1:N2]

    return Y_r
end
