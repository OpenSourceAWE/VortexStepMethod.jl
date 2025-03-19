# replace va_norm_array = norm.(eachrow(solver.sol.va_array)) with a for loop

# testcase
using Test
using LinearAlgebra  # for the norm function

# Define a mock Solver struct
struct MockSolver
    sol::NamedTuple
end

function calc_norm_array!(va_norm_array, va_array)
    for i in 1:size(va_array, 1)
        va_norm_array[i] = norm(view(va_array, i, :))
    end
end

@testset "va_norm_array calculation" begin
    global va_norm_array

    # Create a sample 2D array
    sample_va_array = [
        1.0 2.0 3.0;
        4.0 5.0 6.0;
        7.0 8.0 9.0
    ]

    # Create a mock solver with the sample array
    mock_solver = MockSolver((va_array = sample_va_array,))

    # Calculate va_norm_array
    n = @allocated va_norm_array = norm.(eachrow(mock_solver.sol.va_array))
    println(n)

    va_norm_array2 = zeros(3)
    m = @allocated calc_norm_array!(va_norm_array2, sample_va_array)
    println(m)

    # Expected results (calculated manually)
    expected_norms = [
        sqrt(1^2 + 2^2 + 3^2),
        sqrt(4^2 + 5^2 + 6^2),
        sqrt(7^2 + 8^2 + 9^2)
    ]

    # Test the results
    @test length(va_norm_array) == size(sample_va_array, 1)
    @test va_norm_array ≈ expected_norms atol=1e-10

    # Test individual values
    @test va_norm_array[1] ≈ norm(sample_va_array[1, :]) atol=1e-10
    @test va_norm_array[2] ≈ norm(sample_va_array[2, :]) atol=1e-10
    @test va_norm_array[3] ≈ norm(sample_va_array[3, :]) atol=1e-10

    @test length(va_norm_array2) == size(sample_va_array, 1)
    @test va_norm_array2 ≈ expected_norms atol=1e-10

    # Test individual values
    @test va_norm_array2[1] ≈ norm(sample_va_array[1, :]) atol=1e-10
    @test va_norm_array2[2] ≈ norm(sample_va_array[2, :]) atol=1e-10
    @test va_norm_array2[3] ≈ norm(sample_va_array[3, :]) atol=1e-10
end
nothing