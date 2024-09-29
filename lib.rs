// src/lib.rs

pub mod quantum_state;
pub mod quantum_operations;
pub mod reservoir;
pub mod training;
pub mod hardware;
pub mod visualization;

// Re-export important items for easier use
pub use quantum_state::*;
pub use quantum_operations::*;
pub use reservoir::*;
pub use training::*;
pub use hardware::*;
pub use visualization::*;

/// Ryutan: A Quantum Reservoir Computing framework in Rust
///
/// This library provides tools for implementing and simulating
/// quantum reservoir computing systems, including quantum state
/// representations, operations, reservoir dynamics, and training algorithms.

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
