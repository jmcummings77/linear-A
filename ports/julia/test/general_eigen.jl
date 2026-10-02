function verify_general(a,g)
    n=a.rows
    @test length(g.values)==n
    @test size(g.vectors)==(n,n)
    @test issorted(g.values,by=z->(real(z),imag(z)))
    scale=isempty(a.data) ? 1.0 : maximum(abs,a.data)
    scale==0 && (scale=1.0)
    for j in 1:n
        @test abs(sum(abs2,g.vectors[:,j])-1)<1e-10
        for i in 1:n
            av=sum((a[i,k]/scale)*g.vectors[k,j] for k in 1:n)
            @test abs(av-(g.values[j]/scale)*g.vectors[i,j])<1e-8*max(1,n)
        end
    end
end

@testset "general complex and defective eigenpairs" begin
    for values in (Float64[],[7.0],[0.,-1.,1.,0.],[2.,1.,0.,2.],[0.,1.,0.,0.],zeros(4),[1e300,0.,0.,1e-300],[0.,1e300,-1e-300,0.])
        n=isqrt(length(values));a=Matrix64(n,n,values);before=copy(a.data)
        g=eigen_general(a)
        verify_general(a,g)
        @test a.data==before
        if length(values)==4 && values[2]==1e300
            @test isapprox(g.values[1],-im,atol=1e-12)
            @test isapprox(g.values[2],im,atol=1e-12)
        end
    end
    for scale in (1e300,1e-300,1e-310)
        a=Matrix64(2,2,[0,-scale,scale,0]);g=eigen_general(a)
        verify_general(a,g)
        @test abs(imag(g.values[2])/scale-1)<1e-12
    end
end

@testset "general dense residuals and iteration limits" begin
    for n in 1:12
        a=Matrix64(n,n,[(mod(i*17+n*11,31)-15)/8 for i in 0:n*n-1])
        verify_general(a,eigen_general(a))
    end
    a=Matrix64(4,4,[1,2,3,4,5,6,7,8,2,5,1,3,7,2,8,1])
    @test_throws ErrorException eigen_general(a,max_iterations=1)
    for limit in (0,-1,100001,true,1.5)
        @test_throws ArgumentError eigen_general(a,max_iterations=limit)
    end
    @test_throws DimensionMismatch eigen_general(Matrix64(2,3))
    a.data[1]=NaN
    @test_throws ArgumentError eigen_general(a)
end

@testset "triangular spectra preserve original mixed-scale diagonals" begin
    for values in ([1e300,1.,0.,1e-300],[1e300,0.,1.,1e-300])
        a=Matrix64(2,2,values);g=eigen_general(a)
        @test g.values==ComplexF64[1e-300,1e300]
        verify_general(a,g)
    end
    @test_throws ErrorException eigen_general(Matrix64(2,2,fill(floatmax(Float64),4)))
end
