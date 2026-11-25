use std::num::NonZeroU8;

#[derive(Debug, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
}

#[derive(Debug, Default)]
pub struct AtomTable {
    element: Vec<u8>,
    charge: Vec<i8>,
    n_hydrogens: Vec<u8>,
    isotope: Vec<Option<NonZeroU8>>, // NOTE: keeping u8 for one byte per field, might need u16
    is_aromatic: Vec<bool>,
}

#[derive(Debug, Default)]
pub struct BondTable {
    atom_a: Vec<u32>,
    atom_b: Vec<u32>,
    order: Vec<BondOrder>,
    is_aromatic: Vec<bool>,
}

// pub struct TetrahedralStereo {}
// pub struct DoubelBondStereo {}

#[derive(Debug, Default)]
pub struct LigandGraph {
    atoms: AtomTable,
    bonds: BondTable,

    // NOTE: potential followups if there is time
    // tetrahedral_stereo: Vec<TetrahedralStereo>,
    // double_bond_stereo: Vec<DoubleBondStereo>,
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum LigandGraphError {
    TooManyAtoms,
}

impl LigandGraph {
    pub(crate) fn add_atom(&mut self, element: u8) -> Result<u32, LigandGraphError> {
        let idx = u32::try_from(self.atoms.element.len())
            .map_err(|_| LigandGraphError::TooManyAtoms)?;

        self.atoms.element.push(element);
        self.atoms.charge.push(0);
        self.atoms.n_hydrogens.push(0);
        self.atoms.isotope.push(None);
        self.atoms.is_aromatic.push(false);

        Ok(idx)
    }

    pub(crate) fn add_bond(&mut self, atom_a: u32, atom_b: u32, order: BondOrder) {
        self.bonds.atom_a.push(atom_a);
        self.bonds.atom_b.push(atom_b);
        self.bonds.order.push(order);
        self.bonds.is_aromatic.push(false);
    }

    pub(crate) fn set_hydrogen_count(&mut self, atom: usize, count: u8) {
        self.atoms.n_hydrogens[atom] = count;
    }
}

#[cfg(test)]
mod tests {
    use super::BondOrder::{Double, Single};
    use crate::smiles::smiles_to_ligand_graph;

    #[test]
    fn parse_acetone_in_input_order() {
        let graph = smiles_to_ligand_graph("CC(=O)C")
            .expect("acetone should parse");

        assert_eq!(graph.atoms.element, [6, 6, 8, 6]);
        assert_eq!(graph.atoms.n_hydrogens, [3, 0, 0, 3]);
        assert_eq!(graph.bonds.atom_a, [0, 1, 1]);
        assert_eq!(graph.bonds.atom_b, [1, 2, 3]);
        assert_eq!(graph.bonds.order, [Single, Double, Single]);
    }
}

