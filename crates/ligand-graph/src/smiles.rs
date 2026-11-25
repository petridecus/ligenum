use crate::graph::{BondOrder, LigandGraph, LigandGraphError};

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum SmilesParserError {
    UnsupportedCharacter {
        byte_pos: usize,
        character: char,
    },
    InvalidSyntax {
        byte_pos: usize,
        message: &'static str,
    },
    ValenceExceeded {
        atom_index: u32,
        supported_valence: u8,
        bond_order_sum: usize,
    },
    Graph(LigandGraphError),
}

impl From<LigandGraphError> for SmilesParserError {
    fn from(error: LigandGraphError) -> Self {
        Self::Graph(error)
    }
}

#[derive(Debug)]
struct AtomValence {
    supported: u8,
    bond_order_sum: usize,
}

#[derive(Debug, Default)]
struct SMILESParser {
    graph: LigandGraph,
    current: Option<u32>,
    branches: Vec<u32>,
    pending_order: Option<BondOrder>,
    valences: Vec<AtomValence>,
}

impl SMILESParser {
    fn add_atom(&mut self, ch: char, pos: usize) -> Result<(), SmilesParserError> {
        let (element, supported_valence) = match ch {
            'C' => (6, 4),
            'O' => (8, 2),
            _ => return Err(SmilesParserError::UnsupportedCharacter {byte_pos: pos, character: ch}),
        };

        let child = self.graph.add_atom(element)?;

        self.valences.push(AtomValence {supported: supported_valence, bond_order_sum: 0});

        if let Some(parent) = self.current {
            let order = self.pending_order
                .take()
                .unwrap_or(BondOrder::Single);

            self.connect(parent, child, order);
        }

        self.current = Some(child);
        Ok(())
    }

    fn connect(&mut self, parent: u32, child: u32, order: BondOrder) {
        let weight = match order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        };

        self.graph.add_bond(parent, child, order);
        self.valences[parent as usize].bond_order_sum += weight;
        self.valences[child as usize].bond_order_sum += weight;
    }

    fn open_branch(&mut self, pos: usize, next: Option<char>) -> Result<(), SmilesParserError> {
        let parent = self.current.ok_or(
            SmilesParserError::InvalidSyntax {
                byte_pos: pos,
                message: "attempted to branch without a preceding atom",
            }
        )?;

        if matches!(next, None | Some(')' | '(')) {
            return Err(SmilesParserError::InvalidSyntax {
                byte_pos: pos,
                message: "attempted to create a branch that didn't start with an atom or bond"
            });
        }

        self.branches.push(parent);
        Ok(())
    }

    fn close_branch(&mut self, pos: usize) -> Result<(), SmilesParserError> {
        let parent = self.branches.pop()
            .ok_or(SmilesParserError::InvalidSyntax {
                byte_pos: pos,
                message: "found a closing paren without a matching open",
            })?;

        self.current = Some(parent);
        Ok(())
    }

    fn set_double_bond(&mut self, pos: usize, next: Option<char>) -> Result<(), SmilesParserError> {
        if self.current.is_none() {
            return Err(SmilesParserError::InvalidSyntax {
                byte_pos: pos,
                message: "found '=' without a preceding atom",
            });
        }

        if matches!(next, None | Some('(' | ')' | '=')) {
            return Err(SmilesParserError::InvalidSyntax {
                byte_pos: pos,
                message: "found '=' without a subsequent atom",
            });
        }

        self.pending_order = Some(BondOrder::Double);

        Ok(())
    }

    fn finish(mut self, end_pos: usize) -> Result<LigandGraph, SmilesParserError> {
        if self.current.is_none() {
            return Err(SmilesParserError::InvalidSyntax {
                byte_pos: 0,
                message: "passed in an emtpy SMILES",
            });
        }

        if !self.branches.is_empty() {
            return Err(SmilesParserError::InvalidSyntax {
                byte_pos: end_pos,
                message: "found an open paren without a matching closing",
            });
        }

        for (index, valence) in self.valences.into_iter().enumerate() {
            let atom = index as u32;

            if valence.bond_order_sum > usize::from(valence.supported) {
                return Err(SmilesParserError::ValenceExceeded {
                    atom_index: atom,
                    supported_valence: valence.supported,
                    bond_order_sum: valence.bond_order_sum,
                });
            }

            let hydrogens = valence.supported - valence.bond_order_sum as u8;
            self.graph.set_hydrogen_count(atom as usize, hydrogens);
        }

        Ok(self.graph)
    }
}

pub fn smiles_to_ligand_graph(smiles: &str) -> Result<LigandGraph, SmilesParserError> {
    let mut parser = SMILESParser::default();
    let mut chars = smiles.char_indices().peekable();

    while let Some((pos, ch)) = chars.next() {
        let next = chars.peek().map(|&(_, ch)| ch);

        match ch {
            '(' => parser.open_branch(pos, next)?,
            ')' => parser.close_branch(pos)?,
            '=' => parser.set_double_bond(pos, next)?,
            _ => parser.add_atom(ch, pos)?,
        }
    }

    parser.finish(smiles.len())
}
