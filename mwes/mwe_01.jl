# Replace va_dist = norm.(eachrow(solver.sol._va_dist)) with a for loop

# Testcase that shows that the new function is equivalent to the old, allocating line of code.
using Test
using LinearAlgebra  # for the norm function

# Define a mock Solver struct
struct MockSolver
    sol::NamedTuple
end

function calc_norm_dist!(va_dist, va_vec_dist)
    for i in 1:size(va_vec_dist, 1)
        va_dist[i] = norm(view(va_vec_dist, i, :))
    end
end

@testset "va_dist calculation" begin
    global va_dist

    # Create a sample 2D array
    sample_va_vec_dist = [
        1.0 2.0 3.0;
        4.0 5.0 6.0;
        7.0 8.0 9.0
    ]

    # Create a mock solver with the sample array
    mock_solver = MockSolver((va_vec_dist = sample_va_vec_dist,))

    # Calculate va_dist
    n = @allocated va_dist = norm.(eachrow(mock_solver.sol.va_vec_dist))
    println(n)

    va_dist2 = zeros(3)
    m = @allocated calc_norm_dist!(va_dist2, sample_va_vec_dist)
    println(m)

    # Expected results (calculated manually)
    expected_norms = [
        sqrt(1^2 + 2^2 + 3^2),
        sqrt(4^2 + 5^2 + 6^2),
        sqrt(7^2 + 8^2 + 9^2)
    ]

    # Test the results
    @test length(va_dist) == size(sample_va_vec_dist, 1)
    @test va_dist ≈ expected_norms atol=1e-10

    # Test individual values
    @test va_dist[1] ≈ norm(sample_va_vec_dist[1, :]) atol=1e-10
    @test va_dist[2] ≈ norm(sample_va_vec_dist[2, :]) atol=1e-10
    @test va_dist[3] ≈ norm(sample_va_vec_dist[3, :]) atol=1e-10

    @test length(va_dist2) == size(sample_va_vec_dist, 1)
    @test va_dist2 ≈ expected_norms atol=1e-10

    # Test individual values
    @test va_dist2[1] ≈ norm(sample_va_vec_dist[1, :]) atol=1e-10
    @test va_dist2[2] ≈ norm(sample_va_vec_dist[2, :]) atol=1e-10
    @test va_dist2[3] ≈ norm(sample_va_vec_dist[3, :]) atol=1e-10
end
nothing