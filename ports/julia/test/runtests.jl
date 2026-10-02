using Test
include(joinpath(@__DIR__, "..", "src", "LinearAMatrices.jl"))
using .LinearAMatrices

@testset "owned storage and rectangular arithmetic" begin
    values = Float64[1, 2, 3, 4, 5, 6]
    a = Matrix64(2, 3, values)
    values[1] = 99
    @test rowmajor(transpose(a)) == [1, 4, 2, 5, 3, 6]
    @test rowmajor(a * Matrix64(3, 2, [7, 8, 9, 10, 11, 12])) == [58, 64, 139, 154]
    @test rowmajor(scale(a, 0.5)) == [0.5, 1, 1.5, 2, 2.5, 3]
    @test rowmajor(a + a - a) == rowmajor(a)
    @test rowmajor(identitymatrix(2) * a) == rowmajor(a)
    copied = copy(a)
    copied[1, 1] = 9
    extractedrow, extractedcol, exported = row(a, 1), column(a, 2), rowmajor(a)
    extractedrow[1] = extractedcol[1] = exported[1] = 9
    @test rowmajor(a) == [1, 2, 3, 4, 5, 6]
    @test size(a) == (2, 3)
end

@testset "known determinants and triangular structure" begin
    cases = [(1, [-7], -7), (2, [1, 2, 3, 4], -2), (3, [6, 1, 1, 4, -2, 5, 2, 8, 7], -306),
             (4, [3, 2, 0, 1, 4, 0, 1, 2, 3, 0, 2, 1, 9, 2, 3, 1], 24),
             (3, [0, 1, 0, 0, 0, 1, 1, 0, 0], 1), (2, [1, 2, 2, 4], 0)]
    for (size, values, expected) in cases
        a = Matrix64(size, size, values)
        @test determinant(a) ≈ expected atol=1e-10
        @test rowmajor(a) == values
    end
    a = Matrix64(3, 3, [0.5, 7, -1, 0, -1.5, 4, 0, 0, 2.5])
    @test triangular(a) == (upper=true, lower=false)
    @test triangular(transpose(a)) == (upper=false, lower=true)
    @test trace(a) == 1.5
    @test determinant(a) == -1.875
    @test triangular(Matrix64(2, 3)) == (upper=false, lower=false)
end

@testset "empty shapes and explicit errors" begin
    @test determinant(Matrix64()) == 1
    @test trace(Matrix64()) == 0
    @test triangular(Matrix64()) == (upper=true, lower=true)
    @test rowmajor(Matrix64(2, 0) * Matrix64(0, 3)) == zeros(6)
    @test row(Matrix64(2, 0), 2) == Float64[]
    a = Matrix64(1, 1, [1e308])
    @test_throws ArgumentError scale(a, 2)
    @test_throws ArgumentError a[1, 1] = NaN
    @test_throws BoundsError a[0, 1]
    @test_throws ArgumentError Matrix64(-1, 2)
    @test_throws ArgumentError Matrix64(2, 2, [1])
    @test_throws DimensionMismatch a + Matrix64(2, 2)
    @test_throws DimensionMismatch a * Matrix64(2, 1)
    @test_throws DimensionMismatch trace(Matrix64(2, 3))
    @test_throws DimensionMismatch determinant(Matrix64(2, 3))
    @test a[1, 1] == 1e308
end

@testset "CLI protocol" begin
    runner = joinpath(@__DIR__, "..", "runner.jl")
    command = `$(Base.julia_cmd()) --startup-file=no $runner check multiply 1 2 2 1`
    @test strip(read(pipeline(command, stdin=IOBuffer("2 3 4 5")), String)) == "{\"rows\":1,\"cols\":1,\"values\":[23.0]}"
    output = read(`$(Base.julia_cmd()) --startup-file=no $runner bench scale 2 7 2147483646`, String)
    values = [(mod(i * 17 + 2147483646 * 13, 101) - 50) / 16.0 for i in 0:3]
    expected = (values[1] + values[3] + values[4]) * 1.25 * 7
    @test occursin("\"checksum\":$expected", output)
    @test occursin("\"iterations\":7", output)
    invalid = pipeline(ignorestatus(`$(Base.julia_cmd()) --startup-file=no $runner check trace 1 1`), stdin=IOBuffer("1 2"), stdout=devnull, stderr=devnull)
    @test !success(invalid)
end

@testset "determinant algorithms and scaled products" begin
    for algorithm in (:auto, :lu, :cholesky)
        a = Matrix64(2, 2, [4, 2, 2, 3])
        @test determinant(a, algorithm) ≈ 8
        @test rowmajor(a) == [4, 2, 2, 3]
        @test determinant(Matrix64(), algorithm) == 1
        for diagonal in ([1e200, 1e200, 1e-200, 1e-200], [1e-200, 1e-200, 1e200, 1e200])
            a = Matrix64(4,4)
            for i in 1:4
                a[i,i] = diagonal[i]
            end
            @test determinant(a, algorithm) ≈ 1
        end
    end
    for values in ([1,2,0,1], [1,2,2,1], [1,1,1,1])
        @test_throws ArgumentError determinant(Matrix64(2,2,values), :cholesky)
    end
    @test_throws ArgumentError determinant(identitymatrix(2), :unknown)
    for algorithm in (:auto, :lu)
        @test determinant(Matrix64(2,2,[0,1,1,0]), algorithm) == -1
        @test determinant(Matrix64(2,2,[1,2,2,4]), algorithm) == 0
    end
end

@testset "subnormal elimination ratio" begin
    a = Matrix64(3,3,[1e308,1e308,1e-308,1e-308,2e-308,1e308,0,0,1e-308])
    for algorithm in (:auto, :lu)
        @test determinant(a, algorithm) / 1e-308 ≈ 1
    end
end

@testset "symmetric eigenpairs" begin
    function checkeigen(a, expected)
        original = rowmajor(a)
        result = eigen_symmetric(a)
        @test rowmajor(a) == original
        @test issorted(result.values)
        @test size(result.vectors) == size(a)
        scale = max(maximum(abs, original; init=0.0), 1e-300)
        @test result.values ./ scale ≈ expected ./ scale atol=1e-11
        for col in 1:a.rows, row in 1:a.rows
            av = sum((a[row,k]/scale) * result.vectors[k,col] for k in 1:a.cols)
            @test av ≈ (result.values[col]/scale) * result.vectors[row,col] atol=1e-11
        end
        for i in 1:a.rows, j in 1:a.rows
            dot = sum(result.vectors[row,i]*result.vectors[row,j] for row in 1:a.rows)
            @test dot ≈ (i == j ? 1.0 : 0.0) atol=1e-11
        end
    end
    checkeigen(Matrix64(), Float64[])
    checkeigen(Matrix64(3,3), [0.,0.,0.])
    checkeigen(identitymatrix(3), [1.,1.,1.])
    checkeigen(Matrix64(2,2,[2,1,1,2]), [1,3])
    checkeigen(Matrix64(2,2,[0,-2,-2,0]), [-2,2])
    for scale in (1e-300, 1e300)
        checkeigen(Matrix64(2,2,scale .* [2,1,1,2]), scale .* [1,3])
    end
    @test eigen_symmetric(Matrix64(2,2,[1e-300,0,0,1e300])).values == [1e-300,1e300]
    for size in (3,8,16)
        a = Matrix64(size,size)
        for i in 1:size
            a[i,i] = 2
            if i > 1
                a[i,i-1] = a[i-1,i] = -1
            end
        end
        checkeigen(a, [2 - 2cos(k*pi/(size+1)) for k in 1:size])
        size == 8 && @test_throws ErrorException eigen_symmetric(a; max_sweeps=1)
    end
    @test_throws DimensionMismatch eigen_symmetric(Matrix64(2,3))
    @test_throws ArgumentError eigen_symmetric(Matrix64(2,2,[1,1,0,1]))
    @test_throws ArgumentError eigen_symmetric(identitymatrix(2); tolerance=0)
    @test_throws ArgumentError eigen_symmetric(identitymatrix(2); max_sweeps=0)
end

@testset "cross products and active right-handed rotations" begin
    vector(v) = Matrix64(3,1,v)
    ex,ey,ez = vector([1.,0.,0.]), vector([0.,1.,0.]), vector([0.,0.,1.])
    @test rowmajor(cross(ex,ey)) == [0.,0.,1.]
    @test rowmajor(cross(ey,ex)) == [0.,0.,-1.]
    @test size(cross(transpose(ex),ey)) == (1,3)
    a,b = vector([1.,2.,3.]), vector([-4.,5.,-6.])
    c = cross(a,b)
    @test rowmajor(c) == [-27.,-6.,13.]
    @test sum(a.data .* c.data) == sum(b.data .* c.data) == 0
    @test rowmajor(cross(a,a)) == zeros(3)
    @test rowmajor(a) == [1.,2.,3.]
    @test rowmajor(rotationx(pi/2)*ey) ≈ rowmajor(ez) atol=1e-12
    @test rowmajor(rotationy(pi/2)*ez) ≈ rowmajor(ex) atol=1e-12
    @test rowmajor(rotationz(pi/2)*ex) ≈ rowmajor(ey) atol=1e-12
    @test rowmajor(rotation2d(pi/2)*Matrix64(2,1,[1.,0.])) ≈ [0.,1.] atol=1e-12
    for axis in (vector([1.,2.,3.]),vector([1e300,2e300,3e300]),vector([1e-300,2e-300,3e-300]),vector([nextfloat(0.),0.,0.]))
        before = rowmajor(axis)
        unit = vector(before ./ maximum(abs,before))
        for angle in (0.,0.5,-1.3,pi,1e300)
            r = rotation_axis_angle(axis,angle)
            @test rowmajor(transpose(r)*r) ≈ rowmajor(identitymatrix(3)) atol=1e-12
            @test determinant(r) ≈ 1 atol=1e-12
            @test rowmajor(r*unit) ≈ rowmajor(unit) atol=1e-12
            @test rowmajor(r*rotation_axis_angle(axis,-angle)) ≈ rowmajor(identitymatrix(3)) atol=1e-12
            @test rowmajor(cross(r*a,r*b)) ≈ rowmajor(r*c) atol=1e-12
            @test rowmajor(axis) == before
        end
    end
    for bad in (Matrix64(0,0),Matrix64(3,3),Matrix64(2,1))
        @test_throws DimensionMismatch cross(a,bad)
        @test_throws DimensionMismatch rotation_axis_angle(bad,0.5)
    end
    @test_throws ArgumentError rotation_axis_angle(Matrix64(3,1),0.)
    @test_throws ArgumentError cross(vector([1e308,1e308,0.]), vector([0.,1e308,1e308]))
    for angle in (NaN,Inf,-Inf)
        @test_throws ArgumentError rotation2d(angle)
        @test_throws ArgumentError rotation_axis_angle(ex,angle)
    end
    invalid = vector([1.,2.,3.]); invalid.data[1]=NaN
    @test_throws ArgumentError cross(a,invalid)
    @test_throws ArgumentError rotation_axis_angle(invalid,0.5)
    @test rotationz(1e-200)[2,1] == 1e-200
end

include("general_eigen.jl")
