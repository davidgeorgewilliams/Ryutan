// src/quantum_state/mod.rs

use num_complex::Complex64;
use ndarray::{Array1, Array2};

/// Represents a pure quantum state as a complex vector.
///
/// In quantum mechanics, a pure state represents a quantum system in a definite state.
/// It's described by a state vector (or wave function) in a complex Hilbert space.
/// The `amplitudes` here represent the complex coefficients of the state vector in the computational basis.
///
/// Each element in `amplitudes` represents the complex amplitude of a basis state.
/// For an n-qubit system, there are 2^n basis states.
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Chapter 2.2
/// * "Introduction to Quantum Mechanics" by David J. Griffiths, Chapter 3
#[derive(Clone, Debug, PartialEq)]
pub struct PureState {
    amplitudes: Array1<Complex64>,
}

impl PureState {
    /// Create a new pure state with given amplitudes.
    ///
    /// # Arguments
    ///
    /// * `amplitudes` - A vector of complex amplitudes representing the quantum state.
    pub fn new(amplitudes: Vec<Complex64>) -> Self {
        PureState {
            amplitudes: Array1::from(amplitudes),
        }
    }

    /// Get the number of qubits in the state.
    ///
    /// This calculates the number of qubits based on the length of the amplitude vector.
    /// For n qubits, there are 2^n amplitudes.
    ///
    /// # Returns
    ///
    /// The number of qubits in the state.
    pub fn num_qubits(&self) -> usize {
        (self.amplitudes.len() as f64).log2() as usize
    }

    /// Calculate the probability of measuring a specific basis state.
    ///
    /// In quantum mechanics, the probability of measuring a particular basis state
    /// is given by the absolute square of its amplitude (Born rule).
    ///
    /// # Arguments
    ///
    /// * `basis_state` - The index of the basis state to calculate the probability for.
    ///
    /// # Returns
    ///
    /// The probability of measuring the specified basis state.
    ///
    /// # Additional Reading
    ///
    /// * "Quantum Mechanics: Concepts and Applications" by Nouredine Zettili, Chapter 3.4 (for Born rule)
    pub fn probability(&self, basis_state: usize) -> f64 {
        self.amplitudes[basis_state].norm_sqr()
    }

    /// Get a reference to the amplitudes of the quantum state
    ///
    /// This method provides read-only access to the state amplitudes.
    pub fn get_amplitudes(&self) -> &Array1<Complex64> {
        &self.amplitudes
    }

    /// Get a mutable reference to the amplitudes of the quantum state
    ///
    /// This method provides mutable access to the state amplitudes.
    /// Use with caution, as it allows direct modification of the quantum state.
    pub fn get_amplitudes_mut(&mut self) -> &mut Array1<Complex64> {
        &mut self.amplitudes
    }
}

/// Represents a mixed quantum state as a density matrix.
///
/// A mixed state represents a statistical ensemble of pure states. It's described by a density matrix,
/// which is a positive semi-definite Hermitian matrix with trace 1.
///
/// The density matrix represents the most general quantum state, including both pure and mixed states.
/// The diagonal elements represent probabilities, while off-diagonal elements represent quantum coherences.
///
/// # Additional Reading
///
/// * "Quantum Computation and Quantum Information" by Nielsen & Chuang, Chapter 2.4
/// * "Statistical Mechanics" by R.K. Pathria, Chapter 2 (for a statistical interpretation)
#[derive(Clone, Debug, PartialEq)]
pub struct MixedState {
    density_matrix: Array2<Complex64>,
}

impl MixedState {
    /// Create a new mixed state with given density matrix.
    ///
    /// # Arguments
    ///
    /// * `density_matrix` - A 2D array representing the density matrix of the mixed state.
    pub fn new(density_matrix: Array2<Complex64>) -> Self {
        MixedState {
            density_matrix,
        }
    }

    /// Get the number of qubits in the state.
    ///
    /// Similar to the pure state case, this calculates the number of qubits
    /// based on the dimensions of the density matrix.
    ///
    /// # Returns
    ///
    /// The number of qubits in the state.
    pub fn num_qubits(&self) -> usize {
        (self.density_matrix.nrows() as f64).log2() as usize
    }

    /// Calculate the probability of measuring a specific basis state.
    ///
    /// For a mixed state, the probability of measuring a particular basis state
    /// is given by the corresponding diagonal element of the density matrix.
    ///
    /// # Arguments
    ///
    /// * `basis_state` - The index of the basis state to calculate the probability for.
    ///
    /// # Returns
    ///
    /// The probability of measuring the specified basis state.
    ///
    /// # Additional Reading
    ///
    /// * "Quantum Mechanics: Theory and Applications" by Ajoy Ghatak and S. Lokanathan, Chapter 20 (for density matrix formalism)
    pub fn probability(&self, basis_state: usize) -> f64 {
        self.density_matrix[[basis_state, basis_state]].re
    }

    /// Get a reference to the density matrix of the mixed state
    ///
    /// This method provides read-only access to the density matrix.
    pub fn get_density_matrix(&self) -> &Array2<Complex64> {
        &self.density_matrix
    }

    /// Get a mutable reference to the density matrix of the mixed state
    ///
    /// This method provides mutable access to the density matrix.
    /// Use with caution, as it allows direct modification of the quantum state.
    pub fn get_density_matrix_mut(&mut self) -> &mut Array2<Complex64> {
        &mut self.density_matrix
    }
}

/// Tests for the quantum state implementations.
///
/// These tests demonstrate creating and using both pure and mixed states:
/// * The pure state test creates a state |ψ⟩ = (1/√2)|0⟩ + (i/√2)|1⟩, which is an equal superposition of |0⟩ and |1⟩ with a relative phase.
/// * The mixed state test creates a completely mixed single-qubit state, represented by the density matrix ρ = (1/2)|0⟩⟨0| + (1/2)|1⟩⟨1|.
///
/// Both tests verify the correct number of qubits and probabilities.
///
/// # Additional Reading for Quantum Computing Implementations
///
/// * "Quantum Computing: Where Do We Stand?" by Matthias Troyer (arXiv:1801.04307)
/// * "Quantum Algorithm Implementations for Beginners" by Coles et al. (arXiv:1804.03719)
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pure_state() {
        let state = PureState::new(vec![
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(0.0, 1.0 / 2.0_f64.sqrt()),
        ]);
        assert_eq!(state.num_qubits(), 1);
        assert_eq!(state.probability(0), 0.5);
        assert_eq!(state.probability(1), 0.5);
    }

    #[test]
    fn test_mixed_state() {
        let density_matrix = Array2::from_shape_vec(
            (2, 2),
            vec![
                Complex64::new(0.5, 0.0), Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0), Complex64::new(0.5, 0.0),
            ],
        ).unwrap();
        let state = MixedState::new(density_matrix);
        assert_eq!(state.num_qubits(), 1);
        assert_eq!(state.probability(0), 0.5);
        assert_eq!(state.probability(1), 0.5);
    }
}