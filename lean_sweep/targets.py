"""Bridging-lemma targets for the tactic sweep (CLAUDE_RH_MAP.md section 5).

Each target: name, map (which item in CLAUDE_RH_MAP.md §5), stmt (Lean, after `theorem NAME`),
and `tries`: handwritten tactic scripts, tried in addition to the generic battery.
Statements must be TRUE; the sweep only searches for proofs.
"""

TARGETS = [
    dict(name="fejer_nonneg", map="5.2", stmt=
         "(l γ : ℝ) (hl : 0 ≤ l) : 0 ≤ l * (Real.sin (l * γ / 2) / (l * γ / 2)) ^ 2",
         tries=["positivity"]),

    dict(name="H_closed_form", map="5.1", stmt=
         "(l : ℝ) (z : ℂ) (hl : l ≠ 0) (hz : z ≠ 0) : "
         "2 * (Complex.cosh ((l : ℂ) * z) - 1) / ((l : ℂ) * z ^ 2) = "
         "(l : ℂ) * (Complex.sinh ((l : ℂ) * z / 2) / ((l : ℂ) * z / 2)) ^ 2",
         tries=[
             "have h2 : ((l : ℂ) * z) = 2 * ((l : ℂ) * z / 2) := by ring\n"
             "rw [h2, Complex.cosh_two_mul]\n"
             "have h3 := Complex.cosh_sq_sub_sinh_sq ((l : ℂ) * z / 2)\n"
             "have hl' : (l : ℂ) ≠ 0 := by exact_mod_cast hl\n"
             "field_simp\n"
             "linear_combination (8 * (l:ℂ)) * h3",
             "have hl' : (l : ℂ) ≠ 0 := by exact_mod_cast hl\n"
             "have h3 := Complex.cosh_sq_sub_sinh_sq ((l : ℂ) * z / 2)\n"
             "have h4 : Complex.cosh ((l : ℂ) * z) = Complex.cosh ((l : ℂ) * z / 2) ^ 2 + Complex.sinh ((l : ℂ) * z / 2) ^ 2 := by\n"
             "  rw [← Complex.cosh_two_mul]; ring_nf\n"
             "rw [h4]\n"
             "field_simp\n"
             "linear_combination (-8 * (l:ℂ) * z^2) * h3",
         ]),

    dict(name="H_tail_identity", map="5.4", stmt=
         "(l r : ℝ) (hl : 0 < l) (hr : 0 < r) : "
         "(2 * (Real.cosh (l * r) - 1) / (l * r ^ 2)) * Real.exp (-(l * r)) = "
         "(1 - Real.exp (-(l * r))) ^ 2 / (l * r ^ 2)",
         tries=[
             "rw [Real.cosh_eq]\n"
             "have he : Real.exp (l * r) * Real.exp (-(l * r)) = 1 := by rw [← Real.exp_add]; simp\n"
             "field_simp\n"
             "linear_combination (2:ℝ) * Real.exp (-(l*r)) * he - 2 * he",
             "rw [Real.cosh_eq]\n"
             "have he : Real.exp (l * r) * Real.exp (-(l * r)) = 1 := by rw [← Real.exp_add]; simp\n"
             "field_simp\n"
             "nlinarith [he]",
         ]),

    dict(name="two_box_algebra", map="5.5", stmt=
         "{E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] (f g : E) (h : ‖g‖ = ‖f‖) : "
         "‖f - g‖ ^ 2 = 2 * (‖f‖ ^ 2 - inner ℝ f g)",
         tries=["rw [norm_sub_sq_real, h]; ring", "rw [@norm_sub_sq_real, h]; ring"]),

    dict(name="cos_sum_le", map="5.6", stmt=
         "(n : ℕ) (a γ : Fin n → ℝ) (t : ℝ) (ha : ∀ k, 0 ≤ a k) : "
         "∑ k, a k * Real.cos (γ k * t) ≤ ∑ k, a k",
         tries=["apply Finset.sum_le_sum\nintro k _\nexact mul_le_of_le_one_right (ha k) (Real.cos_le_one _)"]),

    dict(name="fixed_window_constant", map="5.7", stmt=
         "2 * (Real.cosh (Real.log 2 * (1/2)) - 1) / (Real.log 2 * (1/2) ^ 2) = "
         "(6 * Real.sqrt 2 - 8) / Real.log 2",
         tries=[
             "have hs : Real.exp (Real.log 2 * (1/2)) = Real.sqrt 2 := by\n"
             "  rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos (by norm_num)]\n"
             "have hs2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)\n"
             "have hl : Real.log 2 ≠ 0 := by positivity\n"
             "rw [Real.cosh_eq, Real.exp_neg, hs]\n"
             "have hsq : Real.sqrt 2 ≠ 0 := by positivity\n"
             "field_simp\n"
             "nlinarith [hs2]",
         ]),

    dict(name="unit_norm_one", map="5.8a", stmt=
         "(p : ℝ) (hp : 0 ≤ p) : (p + 1 + 2 * Real.sqrt p) * (p + 1 - 2 * Real.sqrt p) = (p - 1) ^ 2",
         tries=["have := Real.sq_sqrt hp\nnlinarith [this]", "ring_nf\nrw [Real.sq_sqrt hp]\nring"]),

    dict(name="unit_primes_235", map="5.8b", stmt=
         "(p : ℕ) (hp : p.Prime) : (p - 1) ∣ 4 ↔ p = 2 ∨ p = 3 ∨ p = 5",
         tries=[
             "constructor\n"
             "· intro h\n"
             "  have h1 : p - 1 ≤ 4 := Nat.le_of_dvd (by norm_num) h\n"
             "  have h2 : 2 ≤ p := hp.two_le\n"
             "  interval_cases p <;> simp_all (config := {decide := true})\n"
             "· rintro (rfl | rfl | rfl) <;> decide",
             "constructor\n"
             "· intro h\n"
             "  have h1 : p - 1 ≤ 4 := Nat.le_of_dvd (by norm_num) h\n"
             "  have h2 : 2 ≤ p := hp.two_le\n"
             "  have h3 : p ≤ 5 := by omega\n"
             "  interval_cases p <;> first | omega | (exfalso; revert hp; decide) | (revert h; decide)\n"
             "· rintro (rfl | rfl | rfl) <;> decide",
         ]),

    dict(name="cosh_half_log2", map="5.7 helper", stmt=
         "Real.cosh (Real.log 2 / 2) = 3 * Real.sqrt 2 / 4",
         tries=[
             "have hs : Real.exp (Real.log 2 / 2) = Real.sqrt 2 := by\n"
             "  rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos (by norm_num)]; ring_nf\n"
             "have hs2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)\n"
             "rw [Real.cosh_eq, Real.exp_neg, hs]\n"
             "have hsq : Real.sqrt 2 ≠ 0 := by positivity\n"
             "field_simp\n"
             "nlinarith [hs2]",
         ]),

    dict(name="completely_additive_primes", map="5.9", stmt=
         "(f : ℕ → ℝ) (hf : ∀ n m, 0 < n → 0 < m → f (n * m) = f n + f m) (n : ℕ) (hn : 0 < n) : "
         "f n = ∑ p ∈ n.primeFactors, (n.factorization p : ℝ) * f p",
         tries=[
             "induction n using Nat.recOnPosPrimePosCoprime with\n"
             "| hp p k hp hk =>\n"
             "  rw [hp.primeFactors_pow (by omega), Finset.sum_singleton, hp.factorization_pow, Finsupp.single_eq_same]\n"
             "  induction k with\n"
             "  | zero => omega\n"
             "  | succ j ih =>\n"
             "    rcases Nat.eq_zero_or_pos j with rfl | hj\n"
             "    · simp\n"
             "    · rw [pow_succ, hf _ _ (by positivity) hp.pos, ih hj (by positivity)]; push_cast; ring\n"
             "| h0 => omega\n"
             "| h1 =>\n"
             "  have := hf 1 1 one_pos one_pos; simp at this ⊢; linarith\n"
             "| h a b ha hb hab iha ihb =>\n"
             "  rw [hf a b (by omega) (by omega), iha (by omega), ihb (by omega),\n"
             "      Nat.Coprime.primeFactors_mul hab, Finset.sum_union hab.disjoint_primeFactors,\n"
             "      Nat.factorization_mul (by omega) (by omega)]\n"
             "  congr 1 <;> apply Finset.sum_congr rfl <;> intro p hp <;> simp [Finsupp.add_apply]\n"
             "  · rw [Nat.factorization_eq_zero_of_not_dvd]; · simp\n"
             "    exact fun hd => (Nat.disjoint_primeFactors hab).ne_of_mem hp ((Nat.mem_primeFactors.mpr ⟨(Nat.prime_of_mem_primeFactors hp), hd, by omega⟩)) rfl\n"
             "  · rw [Nat.factorization_eq_zero_of_not_dvd]; · simp\n"
             "    exact fun hd => (Nat.disjoint_primeFactors hab).ne_of_mem ((Nat.mem_primeFactors.mpr ⟨(Nat.prime_of_mem_primeFactors hp), hd, by omega⟩)) hp rfl",
         ]),

    dict(name="growth_real_positive_coeffs", map="5.10a", stmt=
         "(n : ℕ) (c l : Fin n → ℝ) (M : ℝ) (hc : ∀ k, 0 < c k) "
         "(hM : ∀ s : ℝ, 0 ≤ s → ∑ k, c k * Real.exp (l k * s) ≤ M) : ∀ k, l k ≤ 0",
         tries=[
             "intro k\n"
             "by_contra hpos\n"
             "push_neg at hpos\n"
             "have hle : ∀ s : ℝ, 0 ≤ s → c k * Real.exp (l k * s) ≤ M := by\n"
             "  intro s hs\n"
             "  refine le_trans ?_ (hM s hs)\n"
             "  exact Finset.single_le_sum (f := fun j => c j * Real.exp (l j * s)) (fun j _ => by have := hc j; positivity) (Finset.mem_univ k)\n"
             "set s := (|M| / c k) / l k with hsdef\n"
             "have hs0 : 0 ≤ s := by have := hc k; positivity\n"
             "have h1 := hle s hs0\n"
             "have h2 : l k * s = |M| / c k := by rw [hsdef]; field_simp\n"
             "rw [h2] at h1\n"
             "have h3 := Real.add_one_le_exp (|M| / c k)\n"
             "have hck := hc k\n"
             "have h4 : c k * (|M| / c k + 1) ≤ c k * Real.exp (|M| / c k) := mul_le_mul_of_nonneg_left h3 hck.le\n"
             "have h5 : c k * (|M| / c k + 1) = |M| + c k := by field_simp\n"
             "have := le_abs_self M\n"
             "linarith",
         ]),

    dict(name="residue_weight_identity", map="5.11", stmt=
         "(L t : ℝ) (hL : L ≠ 0) (ht : t ≠ 0) : "
         "L / Real.pi ^ 2 * Real.sin (Real.pi * t) ^ 2 = "
         "t ^ 2 * (L * (Real.sin (L * (2 * Real.pi * t / L) / 2) / (L * (2 * Real.pi * t / L) / 2)) ^ 2)",
         tries=[
             "have hp : Real.pi ≠ 0 := Real.pi_ne_zero\n"
             "have h : L * (2 * Real.pi * t / L) / 2 = Real.pi * t := by field_simp\n"
             "rw [h]\n"
             "field_simp",
             "have hp : Real.pi ≠ 0 := Real.pi_ne_zero\n"
             "have h : L * (2 * Real.pi * t / L) / 2 = Real.pi * t := by field_simp\n"
             "rw [h]\n"
             "field_simp\n"
             "ring",
         ]),

    dict(name="hadamard_psd", map="5.12", stmt=
         "(n : ℕ) (A B : Matrix (Fin n) (Fin n) ℝ) (hA : A.PosSemidef) (hB : B.PosSemidef) : (A.hadamard B).PosSemidef",
         tries=["exact hA.hadamard hB", "exact Matrix.PosSemidef.hadamard hA hB"]),
]

# Generic battery: every target is also tried with each of these.
BATTERY = [
    "simp", "simp_all", "norm_num", "ring", "ring_nf", "field_simp", "field_simp\nring",
    "linarith", "nlinarith", "positivity", "omega", "decide", "tauto", "aesop", "bound",
    "grind", "norm_num [Real.cosh_eq, Complex.cosh, Complex.sinh]", "intros; positivity",
    "intros; nlinarith [sq_nonneg (1:ℝ)]",
]
SEARCH = ["exact?", "apply?"]
