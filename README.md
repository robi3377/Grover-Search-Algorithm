# Grover's Search Algorithm — Qiskit Implementation

A from-scratch implementation of Grover's quantum search algorithm using Qiskit.
Built as a learning project alongside a classical brute-force comparison.

---

## What is Grover's Algorithm?

Suppose you have an unsorted database of N items and you want to find one specific item.
Classically, the only option is to check items one by one — O(N) in the worst case, O(N/2) on average.

Grover's algorithm (1996) solves the same problem in **O(√N) quantum operations**.
For N = 8 that's roughly 2–3 iterations instead of up to 8 checks. The speedup grows with N.

The key insight is amplitude amplification: instead of searching randomly, the algorithm
repeatedly boosts the probability of the correct answer while suppressing all wrong ones.

---

## How It Works

The algorithm operates on a superposition of all N basis states simultaneously.
Each iteration applies two operations back-to-back:

```
|ψ⟩  →  Oracle  →  Diffuser  →  Oracle  →  Diffuser  →  ...  →  Measure
```

After the optimal number of iterations (≈ π/4 · √N), measuring the quantum state
yields the target with high probability.

### 1. Initialization

Apply a Hadamard gate to every qubit to create an equal superposition:

```
|000⟩  →  H⊗n  →  1/√8 (|000⟩ + |001⟩ + ... + |111⟩)
```

Every element in the database now has equal amplitude 1/√N.

### 2. The Oracle

The oracle is a black-box function that "marks" the target state by flipping its phase.
It does not reveal which state is marked — it just applies −1 to it:

```
Oracle|x⟩ = -|x⟩   if x is the target
Oracle|x⟩ =  |x⟩   otherwise
```

This phase flip is invisible on its own (you can't see phase directly), but the diffuser
turns it into a measurable amplitude difference.

**How the oracle circuit works (manual construction):**

1. Flip every qubit where the target bit is 0, so the target maps to |111...1⟩.
2. Apply a multi-controlled Z gate (MCZ) — this gives −1 phase only to |111...1⟩.
3. Flip the same qubits back to restore the original basis.

The MCZ itself is built as: `H · MCX · H` on the last qubit, which converts the
controlled-X into a controlled-Z via the Hadamard identity HXH = Z.

### 3. The Diffuser (Inversion About the Mean)

The diffuser amplifies the marked state by reflecting all amplitudes about their average.
Because the oracle gave the target a negative amplitude, after reflection it ends up
significantly larger than all others.

Circuit: `H X MCZ X H`

- H maps the uniform superposition |s⟩ back to |0...0⟩
- X maps |0...0⟩ to |1...1⟩
- MCZ applies −1 to |1...1⟩ (which represents |s⟩ in this basis)
- X and H undo the mapping

Geometrically, each Oracle+Diffuser pair rotates the state vector by 2θ toward the
target state, where sin(θ) = 1/√N. After k = round(π/4 · √N) rotations, the state
is nearly aligned with the target.

### 4. Why Too Many Iterations Hurt

The rotation overshoots if you apply too many iterations. The probability oscillates
sinusoidally — you can see this clearly in the iteration experiment plot.

---

## Project Structure

```
.
├── main.py          # full implementation
├── circuit.png      # circuit diagram (generated on run)
├── histogram.png    # measurement results bar chart
├── iterations.png   # P(target) vs. iteration count
└── README.md
```

## Setup

```bash
python -m venv venv
venv\Scripts\activate        # Windows
pip install qiskit qiskit-aer matplotlib pylatexenc
python main.py
```

---

## Classical vs. Quantum Comparison

| Metric          | Classical (linear scan) | Grover's Algorithm |
|-----------------|-------------------------|--------------------|
| Search time     | O(N)                    | O(√N)              |
| N = 8 (steps)   | up to 8, avg 4          | 2 iterations       |
| N = 1 million   | up to 1 000 000         | ~785 iterations    |
| Requires        | nothing special         | a quantum computer |

The algorithm is provably optimal — no quantum algorithm can do better than O(√N)
for unstructured search (BBBV theorem, 1994).

---

## References

- L. K. Grover, "A fast quantum mechanical algorithm for database search," STOC 1996.
- Nielsen & Chuang, *Quantum Computation and Quantum Information*, Chapter 6.
- Qiskit documentation: https://docs.quantum.ibm.com
