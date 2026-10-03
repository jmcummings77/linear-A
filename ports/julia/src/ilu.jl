"""Owned ILU(0), without fill, permutations or diagonal shifts. Fields are read-only by contract."""
struct ILU0
    factors::CSRMatrix
    diagonal::Vector{Int}
    function ILU0(a::CSRMatrix)
        a.rows==a.cols || throw(ArgumentError("ILU0 requires square matrix"))
        f=CSRMatrix(a.rows,a.cols,a.offsets,a.indices,a.values);n=f.rows;d=zeros(Int,n)
        for i in 1:n
            d[i]=csr_find(f,i,i-1)
            d[i]<=f.offsets[i+1] && f.indices[d[i]]==i-1 || throw(ArgumentError("ILU0 requires stored diagonal"))
        end
        for i in 1:n
            for p in f.offsets[i]+1:d[i]-1
                j=f.indices[p]+1;f.values[p]/=f.values[d[j]]
                isfinite(f.values[p]) || throw(OverflowError("nonfinite ILU0 factor"))
                for q in d[j]+1:f.offsets[j+1]
                    k=csr_find(f,i,f.indices[q])
                    if k<=f.offsets[i+1] && f.indices[k]==f.indices[q]
                        f.values[k]-=f.values[p]*f.values[q]
                        isfinite(f.values[k]) || throw(OverflowError("nonfinite ILU0 factor"))
                    end
                end
            end
            f.values[d[i]]!=0 || throw(ArgumentError("zero ILU0 pivot"))
        end
        new(f,d)
    end
end
function ilu_apply(f::ILU0,b)
    a,d=f.factors,f.diagonal;n=a.rows
    validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
    a.cols==n && length(d)==n && length(b)==n && all(isfinite,b) || throw(ArgumentError("invalid ILU0 vector or factor"))
    for i in 1:n
        a.offsets[i]<d[i]<=a.offsets[i+1] && a.indices[d[i]]==i-1 && a.values[d[i]]!=0 || throw(ArgumentError("invalid ILU0 diagonal"))
    end
    x=Float64.(b)
    for i in 1:n
        for p in a.offsets[i]+1:d[i]-1;x[i]-=a.values[p]*x[a.indices[p]+1];end
        isfinite(x[i]) || throw(OverflowError("nonfinite ILU0 solve"))
    end
    for i in n:-1:1
        for p in d[i]+1:a.offsets[i+1];x[i]-=a.values[p]*x[a.indices[p]+1];end
        x[i]/=a.values[d[i]];isfinite(x[i]) || throw(OverflowError("nonfinite ILU0 solve"))
    end
    x
end
