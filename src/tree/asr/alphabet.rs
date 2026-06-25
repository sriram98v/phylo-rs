/// Trait defining a biological alphabet for ancestral sequence reconstruction.
pub trait Alphabet: Copy + Clone + Sized + Send + Sync {
    /// Number of states in the alphabet (e.g., 4 for DNA, 20 for protein).
    const N_STATES: usize;
    /// Canonical states in order (e.g., b"ACGT").
    const CANONICAL: &'static [u8];
    /// The gap character.
    const GAP: u8;

    /// Returns a probability profile for a given character.
    /// Returns None if the character is not recognized.
    fn profile(c: u8) -> Option<Vec<f64>>;

    /// Returns the index of the state for a given canonical character.
    fn index_of(c: u8) -> Option<usize>;

    /// Returns the canonical character for a given state index.
    fn char_of(i: usize) -> u8;
}

/// DNA/RNA nucleotide alphabet (4 states: A, C, G, T).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Nucleotide;

impl Alphabet for Nucleotide {
    const N_STATES: usize = 4;
    const CANONICAL: &'static [u8] = b"ACGT";
    const GAP: u8 = b'-';

    fn index_of(c: u8) -> Option<usize> {
        match c {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' | b'U' => Some(3),
            _ => None,
        }
    }

    fn char_of(i: usize) -> u8 {
        Self::CANONICAL[i]
    }

    fn profile(c: u8) -> Option<Vec<f64>> {
        match c {
            b'A' => Some(vec![1.0, 0.0, 0.0, 0.0]),
            b'C' => Some(vec![0.0, 1.0, 0.0, 0.0]),
            b'G' => Some(vec![0.0, 0.0, 1.0, 0.0]),
            b'T' | b'U' => Some(vec![0.0, 0.0, 0.0, 1.0]),
            b'R' => Some(vec![0.5, 0.0, 0.5, 0.0]), // A, G
            b'Y' => Some(vec![0.0, 0.5, 0.0, 0.5]), // C, T
            b'S' => Some(vec![0.0, 0.5, 0.5, 0.0]), // C, G
            b'W' => Some(vec![0.5, 0.0, 0.0, 0.5]), // A, T
            b'K' => Some(vec![0.0, 0.0, 0.5, 0.5]), // G, T
            b'M' => Some(vec![0.5, 0.5, 0.0, 0.0]), // A, C
            b'B' => Some(vec![0.0, 1.0/3.0, 1.0/3.0, 1.0/3.0]), // C, G, T
            b'D' => Some(vec![1.0/3.0, 0.0, 1.0/3.0, 1.0/3.0]), // A, G, T
            b'H' => Some(vec![1.0/3.0, 1.0/3.0, 0.0, 1.0/3.0]), // A, C, T
            b'V' => Some(vec![1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0]), // A, C, G
            b'N' => Some(vec![0.25, 0.25, 0.25, 0.25]),
            Self::GAP => Some(vec![1.0; 4]),
            _ => None,
        }
    }
}

/// Protein amino acid alphabet (20 states).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct AminoAcid;

impl Alphabet for AminoAcid {
    const N_STATES: usize = 20;
    const CANONICAL: &'static [u8] = b"ACDEFGHIKLMNPQRSTVWY";
    const GAP: u8 = b'-';

    fn index_of(c: u8) -> Option<usize> {
        Self::CANONICAL.iter().position(|&x| x == c)
    }

    fn char_of(i: usize) -> u8 {
        Self::CANONICAL[i]
    }

    fn profile(c: u8) -> Option<Vec<f64>> {
        if let Some(idx) = Self::index_of(c) {
            let mut p = vec![0.0; 20];
            p[idx] = 1.0;
            return Some(p);
        }
        match c {
            b'B' => { // Aspartic acid (D) or Asparagine (N)
                let mut p = vec![0.0; 20];
                p[3] = 0.5; // D
                p[11] = 0.5; // N
                Some(p)
            }
            b'Z' => { // Glutamic acid (E) or Glutamine (Q)
                let mut p = vec![0.0; 20];
                p[4] = 0.5; // E
                p[14] = 0.5; // Q
                Some(p)
            }
            b'J' => { // Isoleucine (I) or Leucine (L)
                let mut p = vec![0.0; 20];
                p[8] = 0.5; // I
                p[9] = 0.5; // L
                Some(p)
            }
            b'X' | Self::GAP => Some(vec![1.0; 20]),
            _ => None,
        }
    }
}
