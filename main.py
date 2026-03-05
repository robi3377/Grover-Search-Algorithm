import sys
import io
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram

# Windows terminals often default to cp1252 which can't render Qiskit's
# box-drawing characters — force UTF-8 for stdout
if hasattr(sys.stdout, 'buffer'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# ── Config ──────────────────────────────────────────────────────────────────
N_QUBITS = 3       # 2^3 = 8-element search space
TARGET   = 5       # the marked element we're searching for (|101⟩)
SHOTS    = 4096


# ── Classical brute-force ────────────────────────────────────────────────────
def classical_search(database, target):
    """Linear scan — O(N) worst case."""
    for i, val in enumerate(database):
        if val == target:
            return i, i + 1   # (index, steps taken)
    return -1, len(database)


# ── Oracle ───────────────────────────────────────────────────────────────────
def phase_oracle(n, target):
    """
    Phase oracle for a single marked state.

    Strategy:
      1. Flip every qubit where the target bit is 0 (so the target maps to |111...1⟩).
      2. Apply a multi-controlled-Z — this gives -1 phase only to |111...1⟩.
      3. Undo the flips to restore the original basis.

    Net effect: |target⟩ → -|target⟩,  all other states unchanged.
    """
    qc = QuantumCircuit(n, name='Oracle')

    # bits[i] is the expected bit for qubit i (LSB first, matches Qiskit ordering)
    bits = format(target, f'0{n}b')[::-1]

    # step 1: condition qubits so that the target state becomes |111...1⟩
    for i, b in enumerate(bits):
        if b == '0':
            qc.x(i)

    # step 2: multi-controlled Z (MCZ) via H + MCX + H on the last qubit
    #   MCX flips qubit n-1 only when all control qubits are |1⟩
    #   sandwiching with H turns that X into a Z (phase flip)
    qc.h(n - 1)
    qc.mcx(list(range(n - 1)), n - 1)
    qc.h(n - 1)

    # step 3: undo the bit flips from step 1
    for i, b in enumerate(bits):
        if b == '0':
            qc.x(i)

    return qc


# ── Diffuser ─────────────────────────────────────────────────────────────────
def diffuser(n):
    """
    Grover diffusion operator — 'inversion about the mean'.

    Geometrically: reflects the state vector about the uniform superposition |s⟩.
    Circuit: H X MCZ X H
      - H maps |s⟩ to |0...0⟩
      - X maps |0...0⟩ to |1...1⟩
      - MCZ applies -1 phase to |1...1⟩ (= |0...0⟩ in original basis)
      - X, H undo the mapping
    """
    qc = QuantumCircuit(n, name='Diffuser')
    qc.h(range(n))
    qc.x(range(n))
    # MCZ on all qubits (same trick: H + MCX + H)
    qc.h(n - 1)
    qc.mcx(list(range(n - 1)), n - 1)
    qc.h(n - 1)
    qc.x(range(n))
    qc.h(range(n))
    return qc


# ── Full Grover circuit ───────────────────────────────────────────────────────
def grover_circuit(n, target, iterations):
    qc = QuantumCircuit(n, n)

    # uniform superposition over all 2^n basis states
    qc.h(range(n))
    qc.barrier()

    oracle_gate  = phase_oracle(n, target).to_gate()
    diffuser_gate = diffuser(n).to_gate()

    for _ in range(iterations):
        qc.append(oracle_gate, range(n))
        qc.barrier()
        qc.append(diffuser_gate, range(n))
        qc.barrier()

    qc.measure(range(n), range(n))
    return qc


# ── Simulate ──────────────────────────────────────────────────────────────────
def run_circuit(qc, shots=SHOTS):
    sim = AerSimulator()
    compiled = transpile(qc, sim)
    result = sim.run(compiled, shots=shots).result()
    return result.get_counts()


# ── Iteration experiment ──────────────────────────────────────────────────────
def iteration_experiment(n, target, max_iter=6):
    """Run Grover's for k = 1..max_iter and record P(target) each time."""
    target_key = format(target, f'0{n}b')   # e.g. '101' for target=5
    probs = []

    for k in range(1, max_iter + 1):
        qc     = grover_circuit(n, target, k)
        counts = run_circuit(qc)
        prob   = counts.get(target_key, 0) / SHOTS
        probs.append(prob)
        print(f"  k={k}  P(|{target_key}⟩) = {prob:.3f}")

    return probs


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    N       = 2 ** N_QUBITS
    t_key   = format(TARGET, f'0{N_QUBITS}b')   # Qiskit bitstring for target

    print("=" * 50)
    print("Grover's Search Algorithm")
    print("=" * 50)
    print(f"Search space : {N} elements  (n={N_QUBITS} qubits)")
    print(f"Target       : index {TARGET}  →  |{t_key}⟩")

    # ── Classical comparison ─────────────────────────────────────────────────
    database = list(range(N))
    idx, steps = classical_search(database, TARGET)
    print(f"\nClassical brute-force: found at index {idx} after {steps} step(s)")
    print(f"  Worst case: {N} steps  |  Average: {N // 2} steps")

    # optimal number of Grover iterations ≈ π/4 · sqrt(N)
    optimal_k = round(np.pi / 4 * np.sqrt(N))
    print(f"\nOptimal Grover iterations: {optimal_k}  (≈ π/4 · √{N})")
    print(f"Quantum speedup: O(√N) vs O(N) classical\n")

    # ── Build and draw circuit ───────────────────────────────────────────────
    qc = grover_circuit(N_QUBITS, TARGET, optimal_k)

    # save circuit diagram as image (avoids terminal encoding issues on Windows)
    fig0 = qc.draw(output='mpl', fold=40, style='clifford')
    fig0.savefig(os.path.join('images', 'circuit.png'), dpi=150, bbox_inches='tight')
    print("Circuit diagram saved to images/circuit.png")

    # also print a text version for the terminal
    try:
        print(qc.draw(output='text', fold=150))
    except UnicodeEncodeError:
        print("(text circuit skipped — open circuit.png for the diagram)")

    # ── Run and show results ─────────────────────────────────────────────────
    counts = run_circuit(qc)

    print("\nTop measurement outcomes:")
    for state, cnt in sorted(counts.items(), key=lambda x: -x[1])[:6]:
        bar = '█' * int(cnt / SHOTS * 40)
        marker = '  ← target' if state == t_key else ''
        print(f"  |{state}⟩  {cnt:4d}/{SHOTS}  {cnt/SHOTS*100:5.1f}%  {bar}{marker}")

    # histogram
    fig1 = plot_histogram(
        counts,
        title=f"Grover's Algorithm — Marked state: |{t_key}⟩  ({optimal_k} iteration(s))",
        figsize=(9, 4),
        color='steelblue',
    )
    fig1.savefig(os.path.join('images', 'histogram.png'), dpi=150, bbox_inches='tight')

    # ── Iteration experiment ─────────────────────────────────────────────────
    print("\nIteration experiment — P(target) vs. number of iterations:")
    probs = iteration_experiment(N_QUBITS, TARGET, max_iter=6)

    fig2, ax = plt.subplots(figsize=(7, 4))
    ax.plot(range(1, 7), probs, 'o-', color='steelblue', linewidth=2, markersize=8, label='P(target)')
    ax.axhline(1 / N, color='gray', linestyle='--', linewidth=1.2, label=f'random guess (1/{N})')
    ax.axvline(optimal_k, color='tomato', linestyle=':', linewidth=1.5, label=f'optimal k={optimal_k}')
    ax.set_xlabel('Number of Grover Iterations')
    ax.set_ylabel('P(measuring target state)')
    ax.set_title('Probability of Finding Target vs. Grover Iterations')
    ax.set_ylim(0, 1.05)
    ax.set_xticks(range(1, 7))
    ax.legend()
    ax.grid(alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join('images', 'iterations.png'), dpi=150, bbox_inches='tight')

    print("\nPlots saved to images/: circuit.png, histogram.png, iterations.png")
