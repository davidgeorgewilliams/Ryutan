// src/reservoir/mod.rs

use crate::quantum_state::{PureState};
use crate::quantum_operations::{QuantumOperation, measure};
use rand::Rng;
use num_complex::Complex64;

/// Represents a Quantum Reservoir for Reservoir Computing.
///
/// A Quantum Reservoir is a quantum system that processes input data through its quantum dynamics.
/// It's characterized by a quantum state and a set of quantum operations that define its evolution.
///
/// In the context of Reservoir Computing, the quantum reservoir serves as a high-dimensional,
/// nonlinear dynamical system that can be used for various machine learning tasks.
///
/// # Additional Reading
///
/// * "Quantum Reservoir Computing: A Review" by Fujii & Nakajima (2017)
/// * "Quantum reservoir computing: a renormalization approach" by Ghosh et al. (2019)
pub struct QuantumReservoir {
    state: PureState,
    operations: Vec<Box<dyn QuantumOperation>>,
    num_qubits: usize,
}

impl QuantumReservoir {
    /// Creates a new QuantumReservoir with the given initial state.
    ///
    /// # Arguments
    ///
    /// * `initial_state` - The initial quantum state of the reservoir
    pub fn new(initial_state: PureState) -> Self {
        let num_qubits = initial_state.num_qubits();
        QuantumReservoir {
            state: initial_state,
            operations: Vec::new(),
            num_qubits,
        }
    }

    /// Adds a quantum operation to the reservoir's evolution sequence.
    ///
    /// # Arguments
    ///
    /// * `operation` - The quantum operation to be added
    pub fn add_operation(&mut self, operation: Box<dyn QuantumOperation>) {
        self.operations.push(operation);
    }

    /// Evolves the quantum state of the reservoir by applying all operations in sequence.
    pub fn evolve(&mut self) {
        for operation in &self.operations {
            operation.apply(&mut self.state);
        }
    }

    /// Returns a reference to the current quantum state of the reservoir.
    pub fn get_state(&self) -> &PureState {
        &self.state
    }

    /// Encodes classical input data into the quantum state of the reservoir.
    ///
    /// This method uses amplitude encoding to represent classical data in the quantum state.
    ///
    /// # Arguments
    ///
    /// * `input` - A vector of floating-point numbers representing the input data
    ///
    /// # Additional Reading
    ///
    /// * "Quantum Machine Learning" by Schuld & Petruccione, Chapter 5
    pub fn encode_input(&mut self, input: &[f64]) {
        assert!(input.len() <= 1 << self.num_qubits, "Input size exceeds reservoir capacity");

        let mut amplitudes = vec![Complex64::new(0.0, 0.0); 1 << self.num_qubits];
        let norm = input.iter().map(|&x| x * x).sum::<f64>().sqrt();

        for (i, &value) in input.iter().enumerate() {
            amplitudes[i] = Complex64::new(value / norm, 0.0);
        }

        self.state = PureState::new(amplitudes);
    }

    /// Performs a measurement on all qubits and returns the result as classical data.
    ///
    /// This method collapses the quantum state and can be used as a readout mechanism.
    ///
    /// # Returns
    ///
    /// A vector of boolean values representing the measurement outcomes for each qubit
    pub fn measure_all(&mut self) -> Vec<bool> {
        let mut results = Vec::with_capacity(self.num_qubits);
        for qubit in 0..self.num_qubits {
            results.push(measure(&mut self.state, qubit));
        }
        results
    }

    /// Computes the expectation values of Pauli-Z operators for all qubits.
    ///
    /// This method provides a non-destructive way to extract information from the quantum state.
    ///
    /// # Returns
    ///
    /// A vector of floating-point numbers representing the expectation values
    pub fn compute_expectation_values(&self) -> Vec<f64> {
        let mut expectations = vec![0.0; self.num_qubits];
        for (i, &amp) in self.state.get_amplitudes().iter().enumerate() {
            for qubit in 0..self.num_qubits {
                let sign = if i & (1 << qubit) == 0 { 1.0 } else { -1.0 };
                expectations[qubit] += sign * amp.norm_sqr();
            }
        }
        expectations
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quantum_operations::{PauliX, Hadamard};

    #[test]
    fn test_reservoir_evolution() {
        let initial_state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.add_operation(Box::new(PauliX));
        reservoir.add_operation(Box::new(Hadamard));

        reservoir.evolve();

        let final_state = reservoir.get_state();
        assert!((final_state.get_amplitudes()[0] - Complex64::new(-1.0 / 2.0_f64.sqrt(), 0.0)).norm() < 1e-10);
        assert!((final_state.get_amplitudes()[1] - Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_input_encoding() {
        let initial_state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.encode_input(&[0.6, 0.8]);

        let encoded_state = reservoir.get_state();
        assert!((encoded_state.get_amplitudes()[0] - Complex64::new(0.6, 0.0)).norm() < 1e-10);
        assert!((encoded_state.get_amplitudes()[1] - Complex64::new(0.8, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_measurement() {
        let initial_state = PureState::new(vec![
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)
        ]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        let measurements = reservoir.measure_all();
        assert_eq!(measurements.len(), 1);
    }

    #[test]
    fn test_expectation_values() {
        let initial_state = PureState::new(vec![
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)
        ]);
        let reservoir = QuantumReservoir::new(initial_state);

        let expectations = reservoir.compute_expectation_values();
        assert_eq!(expectations.len(), 1);
        assert!((expectations[0] - 0.0).abs() < 1e-10);
    }
}
