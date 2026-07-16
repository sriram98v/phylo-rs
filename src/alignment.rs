use crate::error::AsrError;
use std::collections::HashMap;

/// A multiple sequence alignment.
pub struct Alignment {
    /// Map from taxon name to sequence bytes.
    pub seqs: HashMap<String, Vec<u8>>,
    /// Width of the alignment in sites.
    pub width: usize,
}

impl Alignment {
    /// Parses a FASTA formatted byte slice.
    pub fn from_fasta_bytes(data: &[u8]) -> Result<Self, AsrError> {
        let mut seqs = HashMap::new();
        let mut current_id = String::new();
        let mut current_seq = Vec::new();

        let lines = data.split(|&b| b == b'\n');
        for line in lines {
            let line = line.trim_ascii();
            if line.is_empty() {
                continue;
            }

            if line[0] == b'>' {
                if !current_id.is_empty() {
                    if seqs.contains_key(&current_id) {
                        return Err(AsrError::InvalidAlignment(format!(
                            "Duplicate taxon ID: {}",
                            current_id
                        )));
                    }
                    seqs.insert(current_id.clone(), current_seq);
                }
                // Header is everything from '>' to first whitespace
                let header = &line[1..];
                let end = header
                    .iter()
                    .position(|&b| b.is_ascii_whitespace())
                    .unwrap_or(header.len());
                let raw_name = String::from_utf8_lossy(&header[..end]);
                if raw_name.is_empty() {
                    return Err(AsrError::InvalidAlignment(
                        "Empty taxon ID in FASTA header".to_string(),
                    ));
                }
                current_id = raw_name.into_owned();
                current_seq = Vec::new();
            } else {
                // Normalize to uppercase
                let normalized: Vec<u8> = line.iter().map(|&b| b.to_ascii_uppercase()).collect();
                current_seq.extend(normalized);
            }
        }

        if !current_id.is_empty() {
            if seqs.contains_key(&current_id) {
                return Err(AsrError::InvalidAlignment(format!(
                    "Duplicate taxon ID: {}",
                    current_id
                )));
            }
            seqs.insert(current_id, current_seq);
        }

        if seqs.is_empty() {
            return Err(AsrError::InvalidAlignment("Empty FASTA file".to_string()));
        }

        let width = seqs.values().next().unwrap().len();
        if !seqs.values().all(|s| s.len() == width) {
            return Err(AsrError::InvalidAlignment(
                "Ragged alignment: sequences have different lengths".to_string(),
            ));
        }

        Ok(Self { seqs, width })
    }

    /// Helper for reading FASTA from a file (std only).
    pub fn from_fasta_file(path: &std::path::Path) -> std::io::Result<Self> {
        let bytes = std::fs::read(path)?;
        Self::from_fasta_bytes(&bytes)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }

    /// Compresses the alignment into unique patterns with multiplicities.
    pub fn compress_columns(&self) -> CompressedColumns {
        let mut pattern_to_idx = HashMap::new();
        let mut patterns = Vec::new();
        let mut multiplicity = Vec::new();
        let mut site_to_pattern = Vec::with_capacity(self.width);

        // Order of leaves for the ASR process
        let leaf_order: Vec<String> = self.seqs.keys().cloned().collect();

        for i in 0..self.width {
            let mut col = Vec::with_capacity(leaf_order.len());
            for name in &leaf_order {
                col.push(self.seqs[name][i]);
            }

            let idx = *pattern_to_idx.entry(col.clone()).or_insert_with(|| {
                let p_idx = patterns.len();
                patterns.push(col);
                multiplicity.push(0);
                p_idx
            });

            multiplicity[idx] += 1;
            site_to_pattern.push(idx);
        }

        CompressedColumns {
            patterns,
            site_to_pattern,
            multiplicity,
            leaf_order,
        }
    }
}

/// Compressed representation of an alignment for performance.
pub struct CompressedColumns {
    /// The unique sequence patterns found in the alignment.
    pub patterns: Vec<Vec<u8>>,
    /// Mapping from original site index to the unique pattern index.
    pub site_to_pattern: Vec<usize>,
    /// The number of times each unique pattern appears in the alignment.
    pub multiplicity: Vec<usize>,
    /// The fixed order of taxa used for the patterns.
    pub leaf_order: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_parse() {
        let data = b">Seq1\nACGT\n>Seq2\nAGGT\n";
        let aln = Alignment::from_fasta_bytes(data).unwrap();
        assert_eq!(aln.width, 4);
        assert_eq!(aln.seqs["Seq1"], b"ACGT");
        assert_eq!(aln.seqs["Seq2"], b"AGGT");
    }

    #[test]
    fn test_compression() {
        let mut seqs = HashMap::new();
        seqs.insert("S1".to_string(), b"AAT".to_vec());
        seqs.insert("S2".to_string(), b"ATT".to_vec());
        let aln = Alignment { seqs, width: 3 };

        let comp = aln.compress_columns();
        // Col 0: A,A (P0)
        // Col 1: A,T (P1)
        // Col 2: T,T (P2)
        assert_eq!(comp.patterns.len(), 3);
        assert_eq!(comp.multiplicity, vec![1, 1, 1]);
    }
}
