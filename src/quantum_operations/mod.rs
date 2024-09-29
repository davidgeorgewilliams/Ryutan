// src/quantum_operations/mod.rs

use crate::quantum_state::{PureState, MixedState};
use ndarray::{Array2, ArrayView2};
use num_complex::Complex64;
use rand::Rng;

/// Trait representing a quantum operation.
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
        let density_matrix = state.get_density_matrix();
        let new_density_matrix = x_matrix.dot(&density_matrix).dot(&x_matrix.t());
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
        let density_matrix = state.get_density_matrix();
        let new_density_matrix = y_matrix.dot(&density_matrix).dot(&y_matrix.t());
        *state.get_density_matrix_mut() = new_density_matrix;
    }
}

impl QuantumOperation for PauliZ {
    fn apply(&self, state: &mut PureState) {
        let n = state.num_qubits();
        for i in (1..state.get_amplitudes().len()).step_by(2) {
            state.get_amplitudes_mut()[i] = -state.get_amplitudes()[i];
        }
    }

    fn apply_to_mixed(&self, state: &mut MixedState) {
        let z_matrix = Array2::from_shape_vec((2, 2), vec![
            Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0),
        ]).unwrap();
        let density_matrix = state.get_density_matrix();
        let new_density_matrix = z_matrix.dot(&density_matrix).dot(&z_matrix.t());
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
        let density_matrix = state.get_density_matrix();
        let new_density_matrix = h_matrix.dot(&density_matrix).dot(&h_matrix.t());
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
    fn test_measurement() {
        let mut state = PureState::new(vec![Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0), Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)]);
        let mut zero_count = 0;
        let trials = 1000;
        for _ in 0..trials {
            let mut test_state = state.clone();
            if !measure(&mut test_state, 0) {
                zero_count += 1;
            }
        }
        let zero_fraction = zero_count as f64 / trials as f64;
        assert!((zero_fraction - 0.5).abs() < 0.05); // Allow for some statistical fluctuation
    }
}
