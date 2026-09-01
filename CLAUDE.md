# GPPDiscovery — Claude operating instructions

This repo (formerly `weil-decay`, renamed 2026-08-24 to reflect its actual scope) is the
standing discovery workbench for numeric/exploratory research across the Golden Physics
Project's shadow framework — see README.md for the framework in brief and the current
active threads. This file holds the *process* rules: who owns what, where a piece of work
belongs, and how it relates to git state so nothing gets orphaned again.

## Ownership (2026-09-01, Daniel, explicit)

**GPPDiscovery is Claude's. GPPDiscovery2 is Codex's.** Each worker sets its own discovery
repo up however it likes — there is no shared convention to negotiate here, and no reason
to mirror the other's layout.

That cuts both ways:
- Read `GPPDiscovery2` freely (standing permission, any session) and port anything useful.
  Attribute a port honestly — say so in the commit and in the Supabase record.
- **Never push to `GPPDiscovery2`.** A `codex/discovery-workbench` branch also exists here
  on GPPDiscovery; leave it alone too. Writes only ever go to Claude's own side —
  Claude's repos and branches, and `public.*` Supabase tables, never `codex.*`.

## First, every turn: the Claude↔Codex channel

`GoldenPhysicsProject/GPP-bridge` is the coordination repo between the two workers. Read
`CONVERSATION.md` there **at the start of every turn** — it is append-only, and Codex
leaves messages in it that will not reach you any other way. The bridge also carries the
migration and admin guides (`docs/MATHLIB-4.33-UPGRADE.md`,
`docs/GITHUB-ADMIN-VIA-POSTGRES.md`); before declaring a GitHub operation impossible,
check whether the bridge's Postgres→GitHub API route already does it.

## Where a piece of work belongs

Three repos, three jobs. Getting this wrong is how work becomes invisible:

| Kind of work | Home |
|---|---|
| Standalone numerics, scans, parameter fits, literature checks | **here** (`discovery/<thread>/`, plus the root `point.py`/`fit.py`/`efit.py` scan harness) |
| Exploratory work that directly supports a Lean thread in flight | `GPPVerify/discovery/<thread>/` — next to the formalization it feeds |
| An actual provable statement | `GPPVerify` Lean sources, via PR |

The second row is a real convention, not an accident: `GPPVerify/discovery/` currently
holds `cutkosky_weil/`, `local_field_shadow/` and `shadow_ope/`, and
`website/GPPVERIFY.md` points sessions at `discovery/cutkosky_weil/notes.md` before they
extend that thread. Don't migrate those here — the pointer would break and the notes would
sit further from the Lean files they explain. Do keep this table honest: if a thread here
grows a Lean counterpart, note the cross-reference in both directions rather than
duplicating the write-up.

## Workflow: discovery here -> formalization in GPPVerify

1. Numeric/exploratory work happens here (`point.py`, `fit.py`, `efit.py`, scan results in
   `results.jsonl`/`E_results.jsonl`, written up in `RESULTS.md`/`E_RESULTS.md`, longer
   threads under `discovery/<thread>/`).
2. The moment a result is solid enough to state as a real theorem (not "the numerics are
   suggestive" — an actual provable statement), it gets formalized in
   `GoldenPhysicsProject/GPPVerify` (see `website/GPPVERIFY.md` for the Lean-side rules:
   no `sorry`, no axiom asserting an open claim, small PRs, CI-green before merge,
   `lean_tasks`/`lean_results`/`formalization_queue` in Supabase).
3. Do this same-session where possible: discover here, prove there, merge there, record in
   Supabase, come back here for the next question. Don't let a promising numeric result
   wait multiple sessions to become a Lean PR — that gap is exactly how stray branches
   happen.
4. **Nothing in this repo is proved.** A result here is evidence, never a theorem. Say so
   in the write-up — a numeric check that passes at every sampled point is not a proof, and
   the ledger discipline the framework's own primary source keeps (proven / argued / open /
   conjectural) applies to everything written here.

## Git conventions

**Here: commit straight to `main`.** There is no CI gate on numerics, so there is no reason
to hold a scan off `main` — and it matches the `workflow_dispatch` auto-commit convention
already in `.github/workflows/` (`scan.yml`, `efit.yml`). Feature branches on this repo
are the exception, not the default; if you make one, close it the same session.

**On GPPVerify: PR-only, and Claude's standing branch is `claude/workbench`.** Codex's is
`codex/lean-workbench` — never push there, never delete it. Cut short-lived
`claude/<thread>` branches off `claude/workbench` when a thread needs its own PR, and
delete them on merge. Full topology table in `website/GPPVERIFY.md`.

**Ref deletion does not work through the git proxy** (`git push origin :refs/heads/…`
returns HTTP 403). Use the bridge's `DELETE /git/refs/heads/<branch, / as %2F>` route.

## The branch-hygiene rule (why this section exists)

On 2026-08-24 an audit of GPPVerify turned up 14 stray branches, some months old, several
with real proved content that had simply never been merged or looked at again — pure loss,
not because the math was wrong but because no session closed the loop. Two were rescued
(PR #122); most were dead. By 2026-09-01 another 22 had accumulated; all were verified
merged-or-dead, their heads recorded in `public.research_notes` (restorable with
`git push origin <sha>:refs/heads/<branch>`), and deleted. That cleanup should not need a
third round.

**Every session that touches this repo or GPPVerify ends one of two ways for every branch
it created:**
1. The branch's content is merged — via PR, CI-green on the actual head SHA, verified —
   before the session ends, or
2. The branch is explicitly closed: either deleted (if superseded/dead), or its state is
   recorded in this repo's `RESULTS.md`/`E_RESULTS.md` (if it's unfinished discovery work
   worth resuming) so a future session finds it there, not by archaeology through
   `git branch -a`.

A branch that just sits there with no PR and no note is the failure mode. Don't create that
failure mode — close every loop you open, same session if at all possible.

## Session protocol

Record to Supabase (project `dunrgpupddbmzffntwph`) before the session ends: `gpp_results`
for project/ops work, `public.research_notes` for discovery findings and infrastructure
decisions, `lean_results`/`formalization_queue` for the Lean side. A failed gate or a
refuted conjecture is a result — write it up honestly rather than logging nothing. Never
write a success row for unverified work.

## Standing pointer

Owner, credentials (Supabase `gpp_vault` — never in this file or in source), and the
Codex division of labor: see `website/CLAUDE.md`. Lean-side rules: `website/GPPVERIFY.md`.
This repo doesn't duplicate them.
