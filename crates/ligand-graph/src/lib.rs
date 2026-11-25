mod graph;
mod smiles;

use pyo3::exceptions::{PyOverflowError, PyValueError};
use pyo3::prelude::*;

use crate::graph::LigandGraph;
use crate::smiles::{smiles_to_ligand_graph, SmilesParserError};

#[pyclass(name = "LigandGraph", module="ligenum")]
struct PyLigandGraph {
    graph: LigandGraph,
}

#[pyfunction]
fn parse_smiles(smiles: &str) -> PyResult<PyLigandGraph> {
    let graph = smiles_to_ligand_graph(smiles).map_err(python_error)?;
    Ok(PyLigandGraph {graph})
}

fn python_error(error: SmilesParserError) -> PyErr {
    match error {
        SmilesParserError::UnsupportedCharacter { byte_pos, character } => {
            PyValueError::new_err(format!(
                "Unsupported SMILES character {character:?} at byte {byte_pos}"
            ))
        }
        SmilesParserError::InvalidSyntax {byte_pos, message} => {
            PyValueError::new_err(format!("{message} at byte {byte_pos}"))
        }
        SmilesParserError::ValenceExceeded {
            atom_index, supported_valence, bond_order_sum,
        } => PyValueError::new_err(format!(
            "Atom {atom_index} has bond-order sum {bond_order_sum}; \
            its supported valence is {supported_valence}"
        )),
        SmilesParserError::Graph(crate::graph::LigandGraphError::TooManyAtoms) => {
            PyOverflowError::new_err("Atom count exceeds index capacity of the graph")
        }
    }
}

#[pymodule]
fn _ligand_graph(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyLigandGraph>()?;
    module.add_function(wrap_pyfunction!(parse_smiles, module)?)?;
    Ok(())
}
