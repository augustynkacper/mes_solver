using Plots

if size(ARGS,1)==0
    N=10
else
    N=parse(Int64, ARGS[1])
end

print("Equation will be solved for N=",N,"...")


# integral from function f, on range (a, b)
function integral(f, a::Float64, b::Float64)::Float64
    div::Int32 = 10
    step::Float64 = (b-a)/div
    res::Float64 = 0.0

    for j in 0:div-step
        res += f(a+(j+0.5)*step)*step
    end 

    return res
end

# returns point x_i
function x_i(i::Int64)::Float64
    return 2*i/N
end


# returns e_i(x)
function e(x, i::Int64)::Float64
    if x > x_i(i-1) && x<=x_i(i)
        return N/2*x - i + 1
    elseif x>x_i(i) && x<x_i(i+1)
        return -N/2*x + i + 1
    else
        return 0.0
    end
end

# returns value of (derative of e_i)(x)
function d_e(x::Float64, i::Int64)::Float64
    if x > x_i(i-1) && x<=x_i(i)
        return N/2
    elseif x>x_i(i) && x<x_i(i+1)
        return -N/2
    else
        return 0
    end
end

# returns function that is 
# a product of deratives of u and v
function du_dv(i::Int64, j::Int64)
    return function(x)
        return d_e(x,i)*d_e(x,j)
    end
end


# returns value of B(e_i,e_j)
function B(i::Int64, j::Int64)::Float64
    a = max(0, x_i(i-1), x_i(j-1))
    b = max(x_i(i+1), x_i(j+1))

    return e(0,i)*e(0,j)-integral(du_dv(i,j), a, b)
end

# returns value of L_i(0)
function L(i::Int64)::Float64
    return 20*e(0,i)
end 

# returns vector X which is solution of B*X=L
function result_func(x::Float64, v)
    res = 0.0
    for i=1:N
        res += v[i]*e(x, i-1)
    end
    return res 
end


function solve()    
    # initialize L and B matrix
    L_matrix = zeros(Float64, N+1, 1)
    B_matrix = zeros(Float64, N+1, N+1)

    # elements on main diagonal have the same value
    B_matrix[2,2]=B(1,1)
    
    # count values in matrix B
    for i=1:N+1
        for j=1:i
            if i==j
                #B_matrix[i,j]=B_matrix[2,2]
                B_matrix[i,j]=B(i-1,j-1)
            else
                # B(e_i, e_j)=B(e_j, e_i)
                B_matrix[i,j]=B(i-1,j-1)
                B_matrix[j,i]=B_matrix[i,j]
            end      
        end
    end

    # count values in matrix L
    for i=1:N
        L_matrix[i,1] = L(i-1)
    end
    
    # u(2)=0 condition
    L_matrix[N+1,1]=0
    for i=1:N
        B_matrix[N+1,i]=0
        B_matrix[i,N+1]=0
    end
    B_matrix[N+1, N+1] = 1

    # solution
    X_m = B_matrix\L_matrix
    
    # x axis
    X_plot = [2*i/N for i=0:N]
    # y axis
    Y_plot = [result_func(X_plot[i+1], X_m) for i=0:N]
   
    plot(X_plot, Y_plot, label="funkcja u")
    xlabel!("x")
    ylabel!("y")  
    savefig("plot.png")
    print("\nPlot of function u is saved in file plot.png.")
end


solve()