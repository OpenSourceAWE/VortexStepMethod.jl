#!/bin/bash
# Integration test for end-user setup.
# Run from the repo root: bash test/test_setup.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PASS=0
FAIL=0

pass() { echo "  PASS: $1"; ((++PASS)); }
fail() { echo "  FAIL: $1"; ((++FAIL)); }

cleanup() {
    [[ -n "${TMPDIR_USER:-}" ]] && rm -rf "$TMPDIR_USER"
}
trap cleanup EXIT

# Use xvfb-run if available (headless CI), otherwise run directly
if command -v xvfb-run &>/dev/null; then
    JULIA="xvfb-run julia"
else
    JULIA="julia"
fi

echo "=== End-user setup ==="
TMPDIR_USER=$(mktemp -d)
echo "  tmpdir: $TMPDIR_USER"
cd "$TMPDIR_USER"

# Simulate: mkdir vsm && cd vsm && julia --project=.
# Then: pkg> add VortexStepMethod (use dev for unreleased)
# Then: julia> using VortexStepMethod; VortexStepMethod.install_examples()
$JULIA --project=. -e '
    using Pkg
    Pkg.develop(path="'"$REPO_ROOT"'")
    using VortexStepMethod
    VortexStepMethod.install_examples()
' 2>&1 && pass "install_examples()" || fail "install_examples()"

# Verify example files were copied
for f in menu.jl rectangular_wing.jl V3_kite.jl pyramid_model.jl \
         ram_air_kite.jl stall_model.jl bench.jl cleanup.jl; do
    [[ -f "examples/$f" ]] && pass "copied examples/$f" || fail "copied examples/$f"
done

# Verify data directories were copied
for d in data/ram_air_kite data/TUDELFT_V3_KITE \
         data/pyramid_model; do
    [[ -d "$d" ]] && pass "copied $d/" || fail "copied $d/"
done

# Verify menu.jl loads GLMakie
if grep -q "using GLMakie" "examples/menu.jl"; then
    pass "menu.jl uses GLMakie"
else
    fail "menu.jl uses GLMakie"
fi

# The copied examples activate their own project (@__DIR__), so ensure that
# environment is fully resolved and linked to this checkout before execution.
$JULIA --project=examples -e '
    using Pkg
    Pkg.develop(path="'"$REPO_ROOT"'")
    Pkg.instantiate()
' 2>&1 && pass "instantiated examples environment" || fail "instantiated examples environment"

# Verify correct packages were added
$JULIA --project=. -e '
    using Pkg
    deps = keys(Pkg.project().dependencies)
    for pkg in ("GLMakie", "LaTeXStrings", "CSV", "DataFrames", "Xfoil")
        @assert pkg ∈ deps "Missing dependency: $pkg"
    end
    @assert !("ControlPlots" ∈ deps) "ControlPlots should not be installed"
' 2>&1 && pass "correct packages installed" || fail "correct packages installed"

# Run all examples via menu.jl --run-all (single Julia process)
echo "  Running all examples..."
$JULIA --project=examples examples/menu.jl -- --run-all \
    2>&1 && pass "run all examples" || fail "run all examples"

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
