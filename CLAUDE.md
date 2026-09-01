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
  Attribute a port honestly — say so in the commit and in `CLAUDE_RESEARCH_NOTES.md`.
- **Never push to `GPPDiscovery2`.** A `codex/discovery-workbench` branch also exists here
  on GPPDiscovery; leave it alone too. Writes only ever go to Claude's own side — Claude's
  repos and branches, and the `CLAUDE_*.md` ledgers. Never edit `CODEX_*.md`.

## Every turn: read the channel first, write to it last

`GoldenPhysicsProject/GPP-bridge` is the coordination repo between the two workers, and
both are under the same standing instruction (Daniel, 2026-09-01):

- **Start of every turn:** `cd /home/user/gpp-bridge && git pull -q origin main && tail -80
  CONVERSATION.md`. Append-only; Codex leaves things there that reach you no other way.
- **End of every turn:** append an entry if anything happened the other side can act on,
  then push. Keep it to the signal.
- **Always:** the full detail goes to `CLAUDE_RESEARCH_NOTES.md` in the same repo. Codex's
  is `CODEX_RESEARCH_NOTES.md`. Read the other's freely; **never edit it.** When a channel
  entry starts running long, put the body in the ledger and have the channel point at it by
  date and heading. Extra `.md` files are welcome for standing reference or
  conjectures/hypotheses — one file per job.

The ledgers moved from Supabase to Markdown in the bridge repo on 2026-09-01 (Daniel): Codex
kept getting blocked from Supabase, so the repo is the surface both sides can always reach.
Supabase still holds the historical record and remains the home for ops (`gpp_results`) and
credentials (`gpp_vault`) — but don't write the same research detail to both, or it drifts.

The bridge also carries the migration and admin guides (`docs/MATHLIB-4.33-UPGRADE.md`,
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
`gpp-bridge/rules/GPPVERIFY.md` points sessions at `discovery/cutkosky_weil/notes.md` before they
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
   `GoldenPhysicsProject/GPPVerify` (see `gpp-bridge/rules/GPPVERIFY.md` for the Lean-side rules:
   no `sorry`, no axiom asserting an open claim, small PRs, CI-green before merge,
   the seven CI gates green before merge).
3. Do this same-session where possible: discover here, prove there, merge there, record in
   `CLAUDE_RESEARCH_NOTES.md`, come back here for the next question. Don't let a promising numeric result
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
delete them on merge. Full topology table in `gpp-bridge/rules/GPPVERIFY.md`.

**Ref deletion does not work through the git proxy** (`git push origin :refs/heads/…`
returns HTTP 403). Use the bridge's `DELETE /git/refs/heads/<branch, / as %2F>` route.

## The branch-hygiene rule (why this section exists)

On 2026-08-24 an audit of GPPVerify turned up 14 stray branches, some months old, several
with real proved content that had simply never been merged or looked at again — pure loss,
not because the math was wrong but because no session closed the loop. Two were rescued
(PR #122); most were dead. By 2026-09-01 another 22 had accumulated; all were verified
merged-or-dead, their heads recorded in `CLAUDE_RESEARCH_NOTES.md` (restorable with
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

Before the session ends, write to the bridge: findings and infrastructure decisions to
`CLAUDE_RESEARCH_NOTES.md`, a refuted route to `CLAUDE_CORRECTIONS.md` (with a *how to catch
this shape next time* line — that is the part worth more than the retraction), and anything
Codex can act on to `CONVERSATION.md`. Move a goal in `CLAUDE_RESEARCH_GOALS.md` when it
lands or dies, so nothing dead sits there looking live.

A failed gate or a refuted conjecture **is** a result — write it up honestly rather than
logging nothing. Never write a success row for unverified work.

Supabase still holds the historical record and remains the home for ops (`gpp_results`) and
credentials (`gpp_vault`), but don't write the same research detail to both.

## Standing pointer

All standing rules moved to `GoldenPhysicsProject/GPP-bridge` on 2026-09-01, so Claude and
Codex read one copy instead of two that drift. Owner, credentials (Supabase `gpp_vault` —
never in this file or in source) and the division of labour: `rules/PROJECT.md`. Lean-side
rules: `rules/GPPVERIFY.md`. Dead routes: `CLAUDE_CORRECTIONS.md` — check it before starting
an RH-positivity thread. Goals: `CLAUDE_RESEARCH_GOALS.md`. This repo doesn't duplicate any
of them.
