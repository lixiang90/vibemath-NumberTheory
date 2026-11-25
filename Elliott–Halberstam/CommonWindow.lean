import Mathlib.Algebra.Order.Field.Power
import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Union
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

/-! Counting simultaneous good windows for finitely many exception sets.
The ceiling of bin j is U/(2B)^j in the logarithmic scale. -/
namespace EHCommonWindow
open scoped Classical

noncomputable def ceiling (U B : ℝ) (j : ℕ) : ℝ := U / (2 * B) ^ j

def Occupied (c : ℕ → ℝ) (exceptions : Finset ℝ) (j : ℕ) : Prop :=
  ∃ a ∈ exceptions, c (j + 1) < a ∧ a ≤ c j

private theorem ceiling_antitone {U B : ℝ} (hU : 0 < U) (hB : 1 ≤ B) :
    Antitone (ceiling U B) := by
  intro i j hij
  have h : 1 ≤ 2 * B := by linarith
  exact div_le_div_of_nonneg_left hU.le (pow_pos (by linarith) _)
    (pow_le_pow_right₀ h hij)

private theorem occupied_count (c : ℕ → ℝ) (hc : Antitone c) (E : Finset ℝ) (H : ℕ) :
    (Finset.filter (fun j => Occupied c E j) (Finset.range H)).card ≤ E.card := by
  classical
  let bad := (Finset.range H).filter (Occupied c E)
  have exists_witness (j : ↥bad) : ∃ a ∈ E, c (j.val + 1) < a ∧ a ≤ c j.val :=
    (Finset.mem_filter.mp j.property).2
  let witness (j : ↥bad) := Classical.choose (exists_witness j)
  have spec (j : ↥bad) : witness j ∈ E ∧ c (j.val + 1) < witness j ∧ witness j ≤ c j.val :=
    Classical.choose_spec (exists_witness j)
  change bad.card ≤ E.card
  refine Finset.card_le_card_of_injective
    (f := fun j : ↥bad => (⟨witness j, (spec j).1⟩ : ↥E)) ?_
  intro i j h
  have heq : witness i = witness j := congrArg Subtype.val h
  have separate (v w : ↥bad) (hlt : v.val < w.val) : witness w < witness v :=
    lt_of_le_of_lt ((spec w).2.2.trans (hc (Nat.succ_le_of_lt hlt))) (spec v).2.1
  apply Subtype.ext
  rcases lt_trichotomy i.val j.val with hij | hij | hij
  · exact ((separate i j hij).ne heq.symm).elim
  · exact hij
  · exact ((separate j i hij).ne heq).elim

private theorem clear_bin {U B : ℝ} {j : ℕ} {E : Finset ℝ} (hB : 1 ≤ B)
    (empty : ¬ Occupied (ceiling U B) E j) {a t : ℝ} (ha : a ∈ E)
    (hlo : ceiling U B j / 2 ≤ t) (hhi : t ≤ ceiling U B j) (hat : a ≤ t) :
    B * a ≤ t := by
  by_contra h
  have hpos : 0 < 2 * B := by linarith
  apply empty
  refine ⟨a, ha, ?_, hat.trans hhi⟩
  have hs : ceiling U B (j + 1) = ceiling U B j / (2 * B) := by
    simp [ceiling, pow_succ, div_div]
  rw [hs, div_lt_iff₀ hpos]
  nlinarith [lt_of_not_ge h]

/-- At most one geometric bin per distinct exception is lost to the whole family. -/
theorem simultaneous_windows (family : Finset (Finset ℝ)) (H : ℕ)
    {U B : ℝ} (hU : 0 < U) (hB : 1 ≤ B) :
    ∃ good : Finset ℕ, good ⊆ Finset.range H ∧
      H - (family.biUnion id).card ≤ good.card ∧
      ∀ j ∈ good, ∀ E ∈ family, ∀ a ∈ E, ∀ t : ℝ,
        ceiling U B j / 2 ≤ t → t ≤ ceiling U B j → a ≤ t → B * a ≤ t := by
  classical
  let F := family.biUnion id
  let bad := (Finset.range H).filter (Occupied (ceiling U B) F)
  have hsub : bad ⊆ Finset.range H := Finset.filter_subset _ _
  refine ⟨Finset.range H \ bad, Finset.sdiff_subset, ?_, ?_⟩
  · rw [Finset.card_sdiff_of_subset hsub, Finset.card_range]
    exact Nat.sub_le_sub_left (occupied_count _ (ceiling_antitone hU hB) F H) _
  · intro j hj E hE a ha t hlo hhi hat
    obtain ⟨hjr, hjb⟩ := Finset.mem_sdiff.mp hj
    have empty : ¬ Occupied (ceiling U B) F j :=
      fun h => hjb (Finset.mem_filter.mpr ⟨hjr, h⟩)
    have haF : a ∈ F := Finset.mem_biUnion.mpr ⟨E, hE, ha⟩
    exact clear_bin hB empty haF hlo hhi hat

/-- The common window may be chosen with w at least U / (2 (2B)^R). -/
theorem exists_common_window (family : Finset (Finset ℝ)) (R : ℕ)
    (hR : (family.biUnion id).card ≤ R) {U B : ℝ} (hU : 0 < U) (hB : 1 ≤ B) :
    ∃ w : ℝ, U / (2 * (2 * B)^R) ≤ w ∧ 2 * w ≤ U ∧
      ∀ E ∈ family, ∀ a ∈ E, ∀ t : ℝ, w ≤ t → t ≤ 2 * w → a ≤ t → B * a ≤ t := by
  classical
  obtain ⟨good, hsub, hsize, hclear⟩ := simultaneous_windows family (R + 1) hU hB
  have hpos : 0 < good.card := by omega
  obtain ⟨j, hj⟩ := Finset.card_pos.mp hpos
  have hjle : j ≤ R := Nat.le_of_lt_succ (Finset.mem_range.mp (hsub hj))
  have hmono := ceiling_antitone hU hB
  have htop : ceiling U B j ≤ U := calc
    ceiling U B j ≤ ceiling U B 0 := hmono (Nat.zero_le j)
    _ = U := by simp [ceiling]
  refine ⟨ceiling U B j / 2, ?_, ?_, ?_⟩
  · calc
      U / (2 * (2 * B)^R) = ceiling U B R / 2 := by simp [ceiling, div_div, mul_comm]
      _ ≤ ceiling U B j / 2 := by linarith [hmono hjle]
  · linarith
  · intro E hE a ha t hlo hhi hat
    exact hclear j hj E hE a ha t hlo (by linarith) hat

end EHCommonWindow

#print axioms EHCommonWindow.simultaneous_windows
#print axioms EHCommonWindow.exists_common_window
