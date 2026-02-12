#!/bin/bash
# Integration test for end-user and developer workflows.
# Run from the repo root: bash test/test_workflows.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PASS=0
FAIL=0

pass() { echo "  PASS: $1"; ((++PASS)); }
fail() { echo "  FAIL: $1"; ((++FAIL)); }

cleanup() {
    [[ -n "${TMPDIR_USER:-}" ]] && rm -rf "$TMPDIR_USER"
    [[ -n "${TMPDIR_DEV:-}" ]]  && rm -rf "$TMPDIR_DEV"
}
trap cleanup EXIT

echo "=== End-user workflow ==="
TMPDIR_USER=$(mktemp -d)
echo "  tmpdir: $TMPDIR_USER"

# Simulate: mkdir vsm && cd vsm && julia --project=.
# Use Pkg.develop instead of Pkg.add (unreleased package)
julia --project="$TMPDIR_USER" -e '
    using Pkg
    Pkg.develop(path="'"$REPO_ROOT"'")
    using VortexStepMethod
    VortexStepMethod.install_examples(false)
' 2>&1 && pass "install_examples(false)" || fail "install_examples(false)"

# Check that example files were copied
cd "$TMPDIR_USER"
for f in menu.jl rectangular_wing.jl V3_kite.jl pyramid_model.jl \
         ram_air_kite.jl stall_model.jl bench.jl cleanup.jl; do
    [[ -f "examples/$f" ]] && pass "copied $f" || fail "copied $f"
done

# Check that copied examples use GLMakie, not ControlPlots
for f in menu.jl bench.jl rectangular_wing.jl V3_kite.jl \
         pyramid_model.jl ram_air_kite.jl stall_model.jl; do
    if grep -q "using GLMakie" "examples/$f" && \
       ! grep -q "using ControlPlots" "examples/$f"; then
        pass "$f uses GLMakie"
    else
        fail "$f uses GLMakie"
    fi
done

# Test install_examples(true) adds the right packages
julia --project="$TMPDIR_USER" -e '
    using Pkg, VortexStepMethod
    VortexStepMethod.install_examples(true)
    deps = keys(Pkg.project().dependencies)
    for pkg in ("GLMakie", "LaTeXStrings", "CSV", "DataFrames", "Xfoil")
        @assert pkg ∈ deps "Missing dependency: $pkg"
    end
    @assert !("ControlPlots" ∈ deps) "ControlPlots should not be installed"
' 2>&1 && pass "install_examples(true) adds correct deps" \
       || fail "install_examples(true) adds correct deps"

cd "$REPO_ROOT"

echo ""
echo "=== Developer workflow ==="
TMPDIR_DEV=$(mktemp -d)
echo "  tmpdir: $TMPDIR_DEV"

# Simulate: git clone ... && cd VortexStepMethod.jl
cp -r "$REPO_ROOT" "$TMPDIR_DEV/VortexStepMethod.jl"
cd "$TMPDIR_DEV/VortexStepMethod.jl"

# Check examples/Project.toml has correct deps
if grep -q "GLMakie" examples/Project.toml && \
   grep -q "VortexStepMethod" examples/Project.toml; then
    pass "examples/Project.toml has GLMakie and VortexStepMethod"
else
    fail "examples/Project.toml has GLMakie and VortexStepMethod"
fi

# Simulate: julia --project=examples, then pkg> dev . && instantiate
julia --project=examples -e '
    using Pkg
    Pkg.develop(path=".")
    Pkg.instantiate()
    deps = keys(Pkg.project().dependencies)
    @assert "GLMakie" ∈ deps "Missing GLMakie"
    @assert "VortexStepMethod" ∈ deps "Missing VortexStepMethod"
' 2>&1 && pass "dev . && instantiate succeeds" \
       || fail "dev . && instantiate succeeds"

# Check examples source files use GLMakie
for f in menu.jl bench.jl rectangular_wing.jl V3_kite.jl \
         pyramid_model.jl ram_air_kite.jl stall_model.jl; do
    if grep -q "using GLMakie" "examples/$f" && \
       ! grep -q "using ControlPlots" "examples/$f"; then
        pass "source $f uses GLMakie"
    else
        fail "source $f uses GLMakie"
    fi
done

cd "$REPO_ROOT"

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
