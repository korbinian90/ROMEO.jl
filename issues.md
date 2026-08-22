# ROMEO.jl - issue list

Working notes from the architecture review of 2026-08-19..21. Uncommitted on
purpose: a plan, not a deliverable. Cross-repository items and the ordering
between them are in `issues-stack.md` (X0..X7). Full evidence:
https://claude.ai/code/artifact/1dcb7a4a-4523-46fc-af23-7c5b0ecc3688

Everything below was measured on this machine unless marked *unverified*.
State refers to branch `claude/julia-repos-architecture-review-tj1zzz` (PR #20).

---

## Done on this branch

- **F17 data race in `unwrap_individual!`** (`src/unwrapping.jl:152`). One
  keyword `Dict` was shared across the `Threads.@threads` loop and written by
  every iteration, so an echo could be unwrapped with another echo's magnitude.
  Five identical calls on 4 threads differed by 6.2831855, exactly 2pi; on 1
  thread, or with no magnitude, 0.0. Reached by `romeo -i` (launchers use
  `--threads=auto`) and by `MriResearchTools.qsm_romeo_B0`. Present in every
  release since 1.2.0 (2024-10-03). Fixed by building `args` per iteration;
  `test/threading.jl` fails against the unfixed code (3 passed, 2 failed).
- **CI could not have caught it.** Julia defaults to one thread. `JULIA_NUM_THREADS: 4`
  added to the runtest step, verified from Pkg's `get_threads_spec()` that it
  forwards to the test worker.
- **P1 ASPIRE citation printed unconditionally as well as conditionally.** Users
  who never ran MCPC-3D-S were told to cite the one method with a patent on it;
  users who did got it twice. Now conditional only, verified by running with and
  without `-B`.
- **`src/provenance.jl`** added: `register_citation!`, `write_provenance`,
  `NOTICES`. Lives here because ROMEO has no dependencies, so every package in
  the family can reach it without constraining anything upstream of itself. Uses
  `Libc.strftime` deliberately, to avoid pulling in Dates.
- **`Project.toml`**: version 1.5.0; stale `ArgParse` removed from `[deps]`
  (it was in both `[deps]` and `[weakdeps]` since the January 2026 refactor;
  verified Pkg never installed it); compat on MriResearchTools deliberately left
  at `"2,3"` so the back-edge of the cycle carries no new constraint.

## Open

### R1. `phase2` read-while-write in the same threaded loop
`unwrap_individual!` passes `phase2 = wrapped[:,:,:,e2]`, a slice another thread
may be unwrapping in place. Left alone deliberately: it did not reproduce over
12 runs with 11 echoes, and fixing it means deciding whether `phase2` should be
the original wrapped phase or the already-unwrapped neighbour. That is a
behaviour change, not a bug fix, and it needs the owner's call.
`src/unwrapping.jl`.

### R2. Extract `RomeoApp` into its own package (F15)
The cycle (ROMEO weak-depends on MriResearchTools for the CLI; MriResearchTools
hard-depends on ROMEO) is defused but not removed: nothing crosses it any more,
verified by running the full suite against the *registered* MriResearchTools
3.5.0 with no source overrides. It still means the dependency-free kernel cannot
be released or audited entirely on its own. The README already describes
extraction as the intended architecture. Do it together with the shared CLI
layer (X1), not as a package created for this one problem. This packaging has
now bitten twice: f01596f "fix circular dependency" (January 2026) and PR #20
going red in 35 seconds when the provenance writer was first put in
MriResearchTools.

### R3. `eval` in the echo-time parser (F19/X2)
`ext/RomeoApp/argparse.jl:175,197,212,222`. `romeo -t '(write("PROOF.txt","executed"); [4,8,12])'`
writes the file. See X2 for why it matters and for the shared fix.

### R4. Untyped keyword plumbing (F10)
Options travel as `keyargs...` into a `Dict{Symbol,Any}` and are probed with
`haskey`, so a misspelled keyword is accepted and discarded rather than
rejected. That is the mechanism that hid F1 for as long as it did: ROMEO accepts
`this_keyword_does_not_exist` without complaint. A typed options struct fixes
the class, and is also what makes dispatch statically resolvable for X6.

### R5. Trimming-relevant, if and when X6 is picked up
- `src/region_handling.jl` holds `merge_regions!` and `correct_regions!`, off by
  default at *runtime*. All 14 verifier errors in the first `--trim=safe`
  attempt traced to this one file: the keyword cannot be known false statically,
  so the branch compiles, dragging in `Statistics.median`, whose error path
  calls `repr`, which pulls in the array-display machinery. Making them a
  compile-time choice took the count to zero.
- The two `Threads.@threads` loops had to be serialised for the binary to run;
  trim does not handle `threading_run` yet. A real cost, not a cleanup.
- Result once both were done: 2.29 MB, 39 ms startup, bit-identical output,
  against 3899 ms JIT and a 160 MB bundle today.

### R6. No Aqua, no JET, no formatter
Family-wide (F14), but ROMEO is the one where a reference implementation most
benefits from it.
