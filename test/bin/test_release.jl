using Test

"""
    dry_run_release(top_header, older_header) -> (exit_code, stderr_text)

Run `bin/release --dry-run` against a throwaway repository at version 5.0.0
whose `CHANGELOG.md` holds a section under each of the two headers, with a stub
`gh` on `PATH` so nothing reaches GitHub.
"""
function dry_run_release(top_header, older_header)
    root = mktempdir()
    repo = mkpath(joinpath(root, "repo"))
    write(joinpath(repo, "Project.toml"), """
        name = "Fixture"
        version = "5.0.0"
        """)
    write(joinpath(repo, "CHANGELOG.md"), """
        # Changelog

        $top_header

        ### Added

        - the newer note

        $older_header

        ### Fixed

        - the older note
        """)
    run(`git -C $repo init --quiet`)
    run(`git -C $repo add --all`)
    run(`git -C $repo -c user.email=fixture@example.com -c user.name=Fixture
         commit --quiet -m "fixture"`)

    stub = mkpath(joinpath(root, "stub"))
    gh_stub = joinpath(stub, "gh")
    write(gh_stub, """
        #!/bin/bash
        case "\$1" in
            repo) echo "OpenSourceAWE/Fixture" ;;
        esac
        """)
    chmod(gh_stub, 0o755)

    script = normpath(joinpath(@__DIR__, "..", "..", "bin", "release"))
    stderr_file = joinpath(root, "stderr.txt")
    env = copy(ENV)
    env["PATH"] = stub * ":" * ENV["PATH"]
    command = setenv(`bash $script --dry-run`, env; dir=repo)
    process = run(pipeline(ignorestatus(command); stdout=devnull, stderr=stderr_file))
    return process.exitcode, read(stderr_file, String)
end

@testset "bin/release" begin
    @testset "release refuses an unversioned top section above the last release" begin
        exit_code, stderr_text =
            dry_run_release("## Unreleased", "## Fixture v5.0.0 2026-09-07")
        @test exit_code != 0
        @test occursin("Version mismatch", stderr_text)
    end

    @testset "release accepts a top section naming the package version" begin
        exit_code, stderr_text =
            dry_run_release("## Fixture v5.0.0 2026-09-07", "## Fixture v4.3.1 2026-08-01")
        @test exit_code == 0
        @test isempty(stderr_text)
    end
end
