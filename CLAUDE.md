# GPPDiscovery — Claude operating instructions

This repo (formerly `weil-decay`, renamed 2026-08-24 to reflect its actual scope) is the
standing discovery workbench for numeric/exploratory research across the Golden Physics
Project's shadow framework — see README.md for the framework in brief and the current
active threads. This file holds the *process* rule that was missing before 2026-08-24:
how work here relates to git state, so nothing gets orphaned again.

## The branch-hygiene rule (why this file exists)

On 2026-08-24 an audit of GPPVerify turned up 14 stray branches, some months old, several
with real proved content that had simply never been merged or looked at again — pure loss,
not because the math was wrong but because no session closed the loop. Two were rescued
(PR #122); most were dead (superseded, or descended from a pre-reset history with no path
back to `main`). That cleanup should not need to happen twice.

**Every session that touches this repo or GPPVerify ends one of two ways for every branch
it created:**
1. The branch's content is merged — via PR, CI-green on the actual head SHA, verified —
   before the session ends, or
2. The branch is explicitly closed: either deleted (if superseded/dead), or its state is
   recorded in this repo's `RESULTS.md`/`E_RESULTS.md` (if it's unfinished discovery work
   worth resuming) so a future session finds it there, not by archaeology through `git
   branch -a`.

A branch that just sits there with no PR and no note is the failure mode. Don't create
that failure mode — close every loop you open, same session if at all possible.

## Workflow: discovery here -> formalization in GPPVerify

1. Numeric/exploratory work happens here (`point.py`, `fit.py`, `efit.py`, scan results in
   `results.jsonl`/`E_results.jsonl`, written up in `RESULTS.md`/`E_RESULTS.md`).
2. The moment a result is solid enough to state as a real theorem (not "the numerics are
   suggestive" — an actual provable statement), it gets formalized directly in
   `GoldenPhysicsProject/GPPVerify` (see that repo's own `CLAUDE.md` for the Lean-side
   rules: no `sorry`, no axiom asserting an open claim, small PRs, CI-green before merge,
   `lean_tasks`/`lean_results`/`formalization_queue` in Supabase).
3. Do this same-session where possible: discover here, prove there, merge there, record in
   Supabase, come back here for the next question. Don't let a promising numeric result
   wait multiple sessions to become a Lean PR — that gap is exactly how stray branches
   happen.
4. Commit scan/discovery results back to this repo's `main` directly (matches the existing
   `workflow_dispatch` auto-commit convention in `.github/workflows/`) rather than parking
   them on a feature branch — there's no CI gate on numerics, so there's no reason to hold
   them off `main`.

## Standing pointer

Owner/credentials/Supabase project details: see `GPPVerify/CLAUDE.md` and
`website/CLAUDE.md` (this repo doesn't duplicate them).
