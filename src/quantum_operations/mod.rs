// src/quantum_operations/mod.rs

use crate::quantum_state::{PureState, MixedState};
use ndarray::Array2;
use num_complex::Complex64;
use rand::Rng;

/// Trait representing a single-qubit quantum operation.
///
/// Quantum operations are transformations applied to quantum states. They can be
/// unitary operations (gates) or non-unitary operations (measurements).
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Chapter 4
pub trait QuantumOperation {
    /// Apply the quantum operation to a pure state.
    fn apply(&self, state: &mut PureState);

    /// Apply the quantum operation to a mixed state.
    fn apply_to_mixed(&self, state: &mut MixedState);
}

/// Trait representing a multi-qubit quantum operation.
///
/// Multi-qubit operations act on specific qubits within a quantum state, enabling
/// entanglement and more complex quantum algorithms.
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Chapter 4.3
pub trait MultiQubitOperation {
    /// Apply the multi-qubit operation to a pure state.
    ///
    /// # Arguments
    ///
    /// * `state` - The quantum state to transform
    /// * `qubits` - The indices of the qubits this operation acts upon
    fn apply_to_qubits(&self, state: &mut PureState, qubits: &[usize]);

    /// Apply the multi-qubit operation to a mixed state.
    ///
    /// # Arguments
    ///
    /// * `state` - The mixed quantum state to transform
    /// * `qubits` - The indices of the qubits this operation acts upon
    fn apply_to_mixed_qubits(&self, state: &mut MixedState, qubits: &[usize]);
}

/// Pauli-X gate (bit flip).
///
/// The Pauli-X gate is equivalent to the NOT gate for classical computers
/// with respect to the standard basis |0⟩, |1⟩.
///
/// Matrix representation:
/// [0 1]
/// [1 0]
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Section 4.2
pub struct PauliX;

/// Pauli-Y gate.
///
/// The Pauli-Y gate rotates the state around the Y-axis of the Bloch sphere by π radians.
///
/// Matrix representation:
/// [0 -i]
/// [i  0]
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Section 4.2
pub struct PauliY;

/// Pauli-Z gate (phase flip).
///
/// The Pauli-Z gate leaves the basis state |0⟩ unchanged and maps |1⟩ to -|1⟩.
///
/// Matrix representation:
/// [1  0]
/// [0 -1]
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Section 4.2
pub struct PauliZ;

/// Hadamard gate.
///
/// The Hadamard gate creates an equal superposition of |0⟩ and |1⟩ when applied to either.
///
/// Matrix representation:
/// [1  1] / √2
/// [1 -1] / √2
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Section 4.2
pub struct Hadamard;

/// Controlled-NOT (CNOT) gate.
///
/// The CNOT gate is a two-qubit operation where the second qubit (target) is flipped
/// if and only if the first qubit (control) is in the |1⟩ state. This gate is fundamental
/// for creating entanglement between qubits.
///
/// Matrix representation (4x4):
/// [1 0 0 0]
/// [0 1 0 0]
/// [0 0 0 1]
/// [0 0 1 0]
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Section 4.3
pub struct CNOT;

/// SWAP gate.
///
/// The SWAP gate exchanges the states of two qubits. It can be decomposed into
/// three CNOT gates but is often treated as a primitive operation.
///
/// Matrix representation (4x4):
/// [1 0 0 0]
/// [0 0 1 0]
/// [0 1 0 0]
/// [0 0 0 1]
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Section 4.3
pub struct SWAP;

/// Controlled-Z (CZ) gate.
///
/// The CZ gate applies a phase flip to the target qubit if the control qubit is |1⟩.
/// Unlike CNOT, CZ is symmetric with respect to its two qubits.
///
/// Matrix representation (4x4):
/// [1 0 0  0]
/// [0 1 0  0]
/// [0 0 1  0]
/// [0 0 0 -1]
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Section 4.3
pub struct ControlledZ;

impl QuantumOperation for PauliX {
    fn apply(&self, state: &mut PureState) {
        let n = state.num_qubits();
        for i in 0..(1 << (n - 1)) {
            let i0 = i * 2;
            let i1 = i0 + 1;
            state.get_amplitudes_mut().swap(i0, i1);
        }
    }

    fn apply_to_mixed(&self, state: &mut MixedState) {
        let x_matrix = Array2::from_shape_vec((2, 2), vec![
            Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
        ]).unwrap();

        let density_matrix = state.get_density_matrix().clone();
        let temp = x_matrix.clone().dot(&density_matrix);
        let x_dagger = x_matrix.t().to_owned();
        let new_density_matrix = temp.dot(&x_dagger);
        *state.get_density_matrix_mut() = new_density_matrix;
    }
}

impl QuantumOperation for PauliY {
    fn apply(&self, state: &mut PureState) {
        let n = state.num_qubits();
        let i = Complex64::new(0.0, 1.0);
        for j in 0..(1 << (n - 1)) {
            let j0 = j * 2;
            let j1 = j0 + 1;
            let temp = state.get_amplitudes()[j0];
            state.get_amplitudes_mut()[j0] = -i * state.get_amplitudes()[j1];
            state.get_amplitudes_mut()[j1] = i * temp;
        }
    }

    fn apply_to_mixed(&self, state: &mut MixedState) {
        let y_matrix = Array2::from_shape_vec((2, 2), vec![
            Complex64::new(0.0, 0.0), Complex64::new(0.0, -1.0),
            Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0),
        ]).unwrap();

        let density_matrix = state.get_density_matrix().clone();
        let temp = y_matrix.clone().dot(&density_matrix);
        let y_dagger = y_matrix.t().mapv(|x| x.conj()).to_owned();
        let new_density_matrix = temp.dot(&y_dagger);
        *state.get_density_matrix_mut() = new_density_matrix;
    }
}

impl QuantumOperation for PauliZ {
    fn apply(&self, state: &mut PureState) {
        let _n = state.num_qubits();
        for i in (1..state.get_amplitudes().len()).step_by(2) {
            state.get_amplitudes_mut()[i] = -state.get_amplitudes()[i];
        }
    }

    fn apply_to_mixed(&self, state: &mut MixedState) {
        let z_matrix = Array2::from_shape_vec((2, 2), vec![
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0),
        ]).unwrap();

        let density_matrix = state.get_density_matrix().clone();
        let temp = z_matrix.clone().dot(&density_matrix);
        let z_dagger = z_matrix.t().to_owned();
        let new_density_matrix = temp.dot(&z_dagger);
        *state.get_density_matrix_mut() = new_density_matrix;
    }
}

impl QuantumOperation for Hadamard {
    fn apply(&self, state: &mut PureState) {
        let n = state.num_qubits();
        let factor = Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0);
        for i in 0..(1 << (n - 1)) {
            let i0 = i * 2;
            let i1 = i0 + 1;
            let temp0 = state.get_amplitudes()[i0];
            let temp1 = state.get_amplitudes()[i1];
            state.get_amplitudes_mut()[i0] = factor * (temp0 + temp1);
            state.get_amplitudes_mut()[i1] = factor * (temp0 - temp1);
        }
    }

    fn apply_to_mixed(&self, state: &mut MixedState) {
        let h_matrix = Array2::from_shape_vec((2, 2), vec![
            Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0), Complex64::new(-1.0, 0.0),
        ]).unwrap() / Complex64::new(2.0_f64.sqrt(), 0.0);

        let density_matrix = state.get_density_matrix().clone();
        let temp = h_matrix.clone().dot(&density_matrix);
        let h_dagger = h_matrix.t().to_owned();
        let new_density_matrix = temp.dot(&h_dagger);
        *state.get_density_matrix_mut() = new_density_matrix;
    }
}

impl MultiQubitOperation for CNOT {
    fn apply_to_qubits(&self, state: &mut PureState, qubits: &[usize]) {
        assert_eq!(qubits.len(), 2, "CNOT requires exactly 2 qubits");
        let control = qubits[0];
        let target = qubits[1];
        let n = state.num_qubits();
        assert!(control < n && target < n, "Qubit indices out of range");
        assert_ne!(control, target, "Control and target must be different qubits");

        let control_mask = 1 << (n - 1 - control);
        let target_mask = 1 << (n - 1 - target);

        for i in 0..(1 << n) {
            if i & control_mask != 0 {
                let j = i ^ target_mask;
                if i < j {
                    state.get_amplitudes_mut().swap(i, j);
                }
            }
        }
    }

    fn apply_to_mixed_qubits(&self, state: &mut MixedState, qubits: &[usize]) {
        assert_eq!(qubits.len(), 2, "CNOT requires exactly 2 qubits");
        let n = state.num_qubits();
        let dim = 1 << n;

        let mut cnot_matrix = Array2::zeros((dim, dim));
        for i in 0..dim {
            cnot_matrix[[i, i]] = Complex64::new(1.0, 0.0);
        }

        let control = qubits[0];
        let target = qubits[1];
        let control_mask = 1 << (n - 1 - control);
        let target_mask = 1 << (n - 1 - target);

        for i in 0..dim {
            if i & control_mask != 0 {
                let j = i ^ target_mask;
                cnot_matrix[[i, i]] = Complex64::new(0.0, 0.0);
                cnot_matrix[[i, j]] = Complex64::new(1.0, 0.0);
            }
        }

        let density_matrix = state.get_density_matrix().clone();
        let temp = cnot_matrix.clone().dot(&density_matrix);
        let cnot_dagger = cnot_matrix.t().to_owned();
        let new_density_matrix = temp.dot(&cnot_dagger);
        *state.get_density_matrix_mut() = new_density_matrix;
    }
}

impl MultiQubitOperation for SWAP {
    fn apply_to_qubits(&self, state: &mut PureState, qubits: &[usize]) {
        assert_eq!(qubits.len(), 2, "SWAP requires exactly 2 qubits");
        let qubit1 = qubits[0];
        let qubit2 = qubits[1];
        let n = state.num_qubits();
        assert!(qubit1 < n && qubit2 < n, "Qubit indices out of range");
        assert_ne!(qubit1, qubit2, "Cannot swap a qubit with itself");

        let mask1 = 1 << (n - 1 - qubit1);
        let mask2 = 1 << (n - 1 - qubit2);

        for i in 0..(1 << n) {
            let bit1 = (i & mask1) != 0;
            let bit2 = (i & mask2) != 0;

            if bit1 != bit2 {
                let j = i ^ mask1 ^ mask2;
                if i < j {
                    state.get_amplitudes_mut().swap(i, j);
                }
            }
        }
    }

    fn apply_to_mixed_qubits(&self, state: &mut MixedState, qubits: &[usize]) {
        assert_eq!(qubits.len(), 2, "SWAP requires exactly 2 qubits");
        let n = state.num_qubits();
        let dim = 1 << n;

        let mut swap_matrix = Array2::zeros((dim, dim));

        let qubit1 = qubits[0];
        let qubit2 = qubits[1];
        let mask1 = 1 << (n - 1 - qubit1);
        let mask2 = 1 << (n - 1 - qubit2);

        for i in 0..dim {
            let bit1 = (i & mask1) != 0;
            let bit2 = (i & mask2) != 0;

            let j = if bit1 != bit2 {
                i ^ mask1 ^ mask2
            } else {
                i
            };

            swap_matrix[[i, j]] = Complex64::new(1.0, 0.0);
        }

        let density_matrix = state.get_density_matrix().clone();
        let temp = swap_matrix.clone().dot(&density_matrix);
        let swap_dagger = swap_matrix.t().to_owned();
        let new_density_matrix = temp.dot(&swap_dagger);
        *state.get_density_matrix_mut() = new_density_matrix;
    }
}

impl MultiQubitOperation for ControlledZ {
    fn apply_to_qubits(&self, state: &mut PureState, qubits: &[usize]) {
        assert_eq!(qubits.len(), 2, "Controlled-Z requires exactly 2 qubits");
        let qubit1 = qubits[0];
        let qubit2 = qubits[1];
        let n = state.num_qubits();
        assert!(qubit1 < n && qubit2 < n, "Qubit indices out of range");

        let mask1 = 1 << (n - 1 - qubit1);
        let mask2 = 1 << (n - 1 - qubit2);

        for i in 0..(1 << n) {
            if (i & mask1 != 0) && (i & mask2 != 0) {
                state.get_amplitudes_mut()[i] = -state.get_amplitudes()[i];
            }
        }
    }

    fn apply_to_mixed_qubits(&self, state: &mut MixedState, qubits: &[usize]) {
        assert_eq!(qubits.len(), 2, "Controlled-Z requires exactly 2 qubits");
        let n = state.num_qubits();
        let dim = 1 << n;

        let mut cz_matrix = Array2::zeros((dim, dim));
        for i in 0..dim {
            cz_matrix[[i, i]] = Complex64::new(1.0, 0.0);
        }

        let qubit1 = qubits[0];
        let qubit2 = qubits[1];
        let mask1 = 1 << (n - 1 - qubit1);
        let mask2 = 1 << (n - 1 - qubit2);

        for i in 0..dim {
            if (i & mask1 != 0) && (i & mask2 != 0) {
                cz_matrix[[i, i]] = Complex64::new(-1.0, 0.0);
            }
        }

        let density_matrix = state.get_density_matrix().clone();
        let temp = cz_matrix.clone().dot(&density_matrix);
        let cz_dagger = cz_matrix.t().to_owned();
        let new_density_matrix = temp.dot(&cz_dagger);
        *state.get_density_matrix_mut() = new_density_matrix;
    }
}

/// Perform a measurement on a specific qubit of a pure quantum state.
///
/// Measurement is a fundamental operation in quantum mechanics that collapses
/// the quantum state to one of the basis states. The outcome is probabilistic,
/// based on the amplitudes of the state.
///
/// # Arguments
///
/// * `state` - The pure quantum state to be measured
/// * `qubit` - The index of the qubit to measure
///
/// # Returns
///
/// A boolean representing the measurement outcome (false for |0⟩, true for |1⟩)
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Section 2.2.3
pub fn measure(state: &mut PureState, qubit: usize) -> bool {
    let n = state.num_qubits();
    assert!(qubit < n, "Qubit index out of range");

    let mut prob_one = 0.0;
    let mask = 1 << (n - 1 - qubit);
    for (i, &amp) in state.get_amplitudes().iter().enumerate() {
        if i & mask != 0 {
            prob_one += amp.norm_sqr();
        }
    }

    let mut rng = rand::thread_rng();
    let outcome = rng.gen::<f64>() < prob_one;

    let factor = 1.0 / (if outcome { prob_one } else { 1.0 - prob_one }).sqrt();
    for (i, amp) in state.get_amplitudes_mut().iter_mut().enumerate() {
        if (i & mask != 0) != outcome {
            *amp = Complex64::new(0.0, 0.0);
        } else {
            *amp *= factor;
        }
    }

    outcome
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pauli_x() {
        let mut state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        PauliX.apply(&mut state);
        assert!((state.get_amplitudes()[0].norm() < 1e-10) && (state.get_amplitudes()[1].norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pauli_y() {
        let mut state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        PauliY.apply(&mut state);
        assert!((state.get_amplitudes()[0].norm() < 1e-10) && (state.get_amplitudes()[1] - Complex64::new(0.0, 1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_pauli_z() {
        let mut state = PureState::new(vec![Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)]);
        PauliZ.apply(&mut state);
        assert!((state.get_amplitudes()[0].norm() < 1e-10) && (state.get_amplitudes()[1] + Complex64::new(1.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_hadamard() {
        let mut state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        Hadamard.apply(&mut state);
        assert!((state.get_amplitudes()[0] - Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)).norm() < 1e-10);
        assert!((state.get_amplitudes()[1] - Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_cnot_gate() {
        let mut state = PureState::new(vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]);

        CNOT.apply_to_qubits(&mut state, &[0, 1]);

        assert!((state.get_amplitudes()[0].norm() < 1e-10));
        assert!((state.get_amplitudes()[1].norm() < 1e-10));
        assert!((state.get_amplitudes()[2].norm() < 1e-10));
        assert!((state.get_amplitudes()[3].norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cnot_creates_entanglement() {
        let mut state = PureState::new(vec![
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(0.0, 0.0),
        ]);

        CNOT.apply_to_qubits(&mut state, &[0, 1]);

        let expected_00 = 1.0 / 2.0_f64.sqrt();
        let expected_11 = 1.0 / 2.0_f64.sqrt();

        assert!((state.get_amplitudes()[0].norm() - expected_00).abs() < 1e-10);
        assert!((state.get_amplitudes()[1].norm() < 1e-10));
        assert!((state.get_amplitudes()[2].norm() < 1e-10));
        assert!((state.get_amplitudes()[3].norm() - expected_11).abs() < 1e-10);
    }

    #[test]
    fn test_swap_gate() {
        let mut state = PureState::new(vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]);

        SWAP.apply_to_qubits(&mut state, &[0, 1]);

        assert!((state.get_amplitudes()[0].norm() < 1e-10));
        assert!((state.get_amplitudes()[1].norm() < 1e-10));
        assert!((state.get_amplitudes()[2].norm() - 1.0).abs() < 1e-10);
        assert!((state.get_amplitudes()[3].norm() < 1e-10));
    }

    #[test]
    fn test_controlled_z_gate() {
        let mut state = PureState::new(vec![
            Complex64::new(0.5, 0.0),
            Complex64::new(0.5, 0.0),
            Complex64::new(0.5, 0.0),
            Complex64::new(0.5, 0.0),
        ]);

        ControlledZ.apply_to_qubits(&mut state, &[0, 1]);

        assert!((state.get_amplitudes()[0] - Complex64::new(0.5, 0.0)).norm() < 1e-10);
        assert!((state.get_amplitudes()[1] - Complex64::new(0.5, 0.0)).norm() < 1e-10);
        assert!((state.get_amplitudes()[2] - Complex64::new(0.5, 0.0)).norm() < 1e-10);
        assert!((state.get_amplitudes()[3] - Complex64::new(-0.5, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_measurement() {
        let state = PureState::new(vec![Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0), Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)]);
        let mut zero_count = 0;
        let trials = 1000;
        for _ in 0..trials {
            let mut test_state = state.clone();
            if !measure(&mut test_state, 0) {
                zero_count += 1;
            }
        }
        let zero_fraction = zero_count as f64 / trials as f64;
        assert!((zero_fraction - 0.5).abs() < 0.05);
    }
}