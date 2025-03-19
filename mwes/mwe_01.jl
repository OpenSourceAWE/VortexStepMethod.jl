# replace va_norm_array = norm.(eachrow(solver.sol.va_array)) with a for loop

# testcase
using Test
using LinearAlgebra  # for the norm function

# Define a mock Solver struct
struct MockSolver
    sol::NamedTuple
end

@testset "va_norm_array calculation" begin
    # Create a sample 2D array
    sample_va_array = [
        1.0 2.0 3.0;
        4.0 5.0 6.0;
        7.0 8.0 9.0
    ]

    # Create a mock solver with the sample array
    mock_solver = MockSolver((va_array = sample_va_array,))

    # Calculate va_norm_array
    va_norm_array = norm.(eachrow(mock_solver.sol.va_array))

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
end
nothing