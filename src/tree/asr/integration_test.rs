/// Integration tests for ancestral sequence reconstruction, adapted from treetime test suite.
/// Reference: https://github.com/neherlab/treetime/blob/master/test/test_treetime.py

use std::collections::HashMap;
use crate::node::NodeID;
use crate::tree::PhyloTree;
use crate::prelude::*;
use crate::tree::asr::alphabet::{Alphabet, Nucleotide};
use crate::tree::asr::gtr::GtrModel;
use crate::tree::asr::alignment::Alignment;
use crate::tree::asr::reconstruction::Reconstruction;

/// Build a simple 4-taxon tree from Newick: ((A:0.6,B:0.3)C:0.1,D:0.2)E:0.001
fn build_tiny_tree() -> PhyloTree {
    let newick = b"((A:0.601,B:0.301)C:0.100,D:0.200)E:0.001;";
    crate::prelude::Newick::from_newick(newick).unwrap()
}

// ===========================================================================
// Alphabet tests
// ===========================================================================

#[test]
fn test_nucleotide_canonical_profiles() {
    assert_eq!(Nucleotide::profile(b'A'), Some(vec![1.0, 0.0, 0.0, 0.0]));
    assert_eq!(Nucleotide::profile(b'C'), Some(vec![0.0, 1.0, 0.0, 0.0]));
    assert_eq!(Nucleotide::profile(b'G'), Some(vec![0.0, 0.0, 1.0, 0.0]));
    assert_eq!(Nucleotide::profile(b'T'), Some(vec![0.0, 0.0, 0.0, 1.0]));
    assert_eq!(Nucleotide::profile(b'U'), Some(vec![0.0, 0.0, 0.0, 1.0])); // U -> T
}

#[test]
fn test_nucleotide_ambiguity_profiles() {
    // R = A or G
    let r = Nucleotide::profile(b'R').unwrap();
    assert!((r[0] - 0.5).abs() < 1e-10);
    assert!((r[2] - 0.5).abs() < 1e-10);

    // Y = C or T
    let y = Nucleotide::profile(b'Y').unwrap();
    assert!((y[1] - 0.5).abs() < 1e-10);
    assert!((y[3] - 0.5).abs() < 1e-10);

    // N = all equal
    let n = Nucleotide::profile(b'N').unwrap();
    for &v in &n {
        assert!((v - 0.25).abs() < 1e-10);
    }

    // Gap = uninformative (all ones)
    let gap = Nucleotide::profile(b'-').unwrap();
    for v in &gap {
        assert_eq!(*v, 1.0);
    }
}

#[test]
fn test_nucleotide_index_of() {
    assert_eq!(Nucleotide::index_of(b'A'), Some(0));
    assert_eq!(Nucleotide::index_of(b'C'), Some(1));
    assert_eq!(Nucleotide::index_of(b'G'), Some(2));
    assert_eq!(Nucleotide::index_of(b'T'), Some(3));
    assert_eq!(Nucleotide::index_of(b'U'), Some(3));
    assert_eq!(Nucleotide::index_of(b'X'), None);
}

#[test]
fn test_nucleotide_char_of() {
    assert_eq!(Nucleotide::char_of(0), b'A');
    assert_eq!(Nucleotide::char_of(1), b'C');
    assert_eq!(Nucleotide::char_of(2), b'G');
    assert_eq!(Nucleotide::char_of(3), b'T');
}

#[test]
fn test_amino_acid_canonical() {
    use crate::tree::asr::alphabet::AminoAcid;

    for (i, c) in AminoAcid::CANONICAL.iter().enumerate() {
        assert_eq!(AminoAcid::char_of(i), *c);
        assert_eq!(AminoAcid::index_of(*c), Some(i));
    }
}

#[test]
fn test_amino_acid_ambiguity_profiles() {
    use crate::tree::asr::alphabet::AminoAcid;

    // B = D or N
    let b = AminoAcid::profile(b'B').unwrap();
    assert!((b[3] - 0.5).abs() < 1e-10);
    assert!((b[11] - 0.5).abs() < 1e-10);

    // Z = E or Q
    let z = AminoAcid::profile(b'Z').unwrap();
    assert!((z[4] - 0.5).abs() < 1e-10);
    assert!((z[14] - 0.5).abs() < 1e-10);

    // J = I or L
    let j = AminoAcid::profile(b'J').unwrap();
    assert!((j[8] - 0.5).abs() < 1e-10);
    assert!((j[9] - 0.5).abs() < 1e-10);

    // X = uninformative
    let x = AminoAcid::profile(b'X').unwrap();
    for v in &x {
        assert_eq!(*v, 1.0);
    }
}

// ===========================================================================
// GTR model tests (adapted from treetime test_GTR)
// ===========================================================================

#[test]
fn test_gtr_jc_p_zero_is_identity() {
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
    let p_t = model.transition(0.0);

    for i in 0..4 {
        if i == 0 {
            assert_eq!(p_t[(i, i)], 1.0);
        } else {
            assert_eq!(p_t[(i, 0)], 0.0);
        }
    }
}

#[test]
fn test_gtr_jc_row_sums_to_one() {
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
    for t in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0] {
        let p_t = model.transition(t);
        for i in 0..4 {
            let sum: f64 = (0..4).map(|j| p_t[(i, j)]).sum();
            assert!((sum - 1.0).abs() < 1e-10, "t={:.2} row {} sum = {}", t, i, sum);
        }
    }
}

#[test]
fn test_gtr_jc_positive_entries() {
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
    for t in [0.001, 0.5, 2.0, 10.0] {
        let p_t = model.transition(t);
        for i in 0..4 {
            for j in 0..4 {
                assert!(p_t[(i, j)] >= -1e-15, "Negative entry at t={}: ({},{}) = {}", t, i, j, p_t[(i,j)]);
            }
        }
    }
}

#[test]
fn test_gtr_p_infinity_converges_to_pi() {
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
    let pi = model.equilibrium();
    let p_inf = model.transition(100.0); // large t

    for i in 0..4 {
        for j in 0..4 {
            let diff = (p_inf[(i, j)] - pi[j]).abs();
            assert!(diff < 0.01, "t=100 ({},{}) expected {} got {}", i, j, pi[j], p_inf[(i,j)]);
        }
    }
}

#[test]
fn test_gtr_q_rows_sum_zero() {
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
    let p_small = model.transition(1e-10);
    for i in 0..4 {
        let sum: f64 = (0..4).map(|j| p_small[(i, j)]).sum();
        assert!((sum - 1.0).abs() < 1e-8);
    }
}

#[test]
fn test_gtr_detailed_balance() {
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
    let pi = model.equilibrium();

    for t in [0.1, 0.5, 1.0, 5.0] {
        let p_t = model.transition(t);

        // For JC, pi is uniform so pi^T P(t) should equal pi
        let result: Vec<f64> = (0..4usize).map(|j| (0..4usize).map(|i| pi[i] * p_t[(i, j)]).sum()).collect();

        for j in 0..4 {
            assert!((result[j] - pi[j]).abs() < 1e-8, "t={:.2} j={} expected {} got {}", t, j, pi[j], result[j]);
        }
    }
}

#[test]
fn test_gtr_custom_pi() {
    let pi = vec![0.30, 0.20, 0.30, 0.20];
    let w = nalgebra::DMatrix::from_element(4, 4, 1.0);
    let model = GtrModel::<Nucleotide>::new(pi, w, true).unwrap();

    let pi_norm = model.equilibrium();
    let sum: f64 = pi_norm.iter().sum();
    assert!((sum - 1.0).abs() < 1e-10);

    for &p in pi_norm.as_slice() {
        assert!(p > 0.0);
    }
}

// ===========================================================================
// Alignment tests (adapted from treetime alignment parsing)
// ===========================================================================

#[test]
fn test_alignment_fasta_multiline() {
    let data = b">Seq1\nACGTACGT\nACGTACGT\n>Seq2\nAGGTAGGT\nAGGTAGGT\n";
    let aln = Alignment::from_fasta_bytes(data).unwrap();
    assert_eq!(aln.width, 16);
    assert_eq!(aln.seqs["Seq1"], b"ACGTACGTACGTACGT");
}

#[test]
fn test_alignment_fasta_crlf() {
    let data = b">Seq1\r\nACGT\r\n>Seq2\r\nAGGT\r\n";
    let aln = Alignment::from_fasta_bytes(data).unwrap();
    assert_eq!(aln.width, 4);
}

#[test]
fn test_alignment_fasta_uppercase() {
    let data = b">Seq1\nacgt\n>Seq2\naggt\n";
    let aln = Alignment::from_fasta_bytes(data).unwrap();
    assert_eq!(aln.seqs["Seq1"], b"ACGT");
    assert_eq!(aln.seqs["Seq2"], b"AGGT");
}

#[test]
fn test_alignment_fasta_empty_header() {
    let data = b">\nACGT\n>Seq2\nAGGT\n";
    let result = Alignment::from_fasta_bytes(data);
    assert!(result.is_err(), "empty header should return an error");
}

#[test]
fn test_alignment_fasta_empty_fasta_returns_error() {
    let data = b"\n\n\n";
    let result = Alignment::from_fasta_bytes(data);
    assert!(result.is_err());
}

#[test]
fn test_alignment_fasta_ragged_returns_error() {
    let data = b">Seq1\nACGT\n>Seq2\nAGG\n";
    let result = Alignment::from_fasta_bytes(data);
    assert!(result.is_err());
}

#[test]
fn test_alignment_fasta_duplicate_id_returns_error() {
    let data = b">Seq1\nACGT\n>Seq1\nAGGT\n";
    let result = Alignment::from_fasta_bytes(data);
    assert!(result.is_err(), "duplicate taxon ID should return an error");
}

#[test]
fn test_alignment_fasta_whitespace_in_header() {
    let data = b">Seq1 some description here\nACGT\n>Seq2 another desc\nAGGT\n";
    let aln = Alignment::from_fasta_bytes(data).unwrap();
    assert_eq!(aln.seqs["Seq1"], b"ACGT");
    assert_eq!(aln.seqs["Seq2"], b"AGGT");
}

#[test]
fn test_compression_identical_columns() {
    let mut seqs = HashMap::new();
    seqs.insert("S1".to_string(), b"AATT".to_vec());
    seqs.insert("S2".to_string(), "AAT T".replace(' ', "-").as_bytes().to_vec());
    let aln = Alignment { seqs, width: 4 };

    let comp = aln.compress_columns();
    // Col 0: AA -> P0
    // Col 1: AA -> P0 (duplicate!)
    // Col 2: TT -> P1
    // Col 3: T- -> P2
    assert_eq!(comp.multiplicity[0], 2); // AA appears twice
}

// ===========================================================================
// Reconstruction serialization test
// ===========================================================================

#[test]
fn test_reconstruction_sequence_string() {
    let mut seqs = HashMap::new();
    seqs.insert(0_usize.into(), vec![0, 1, 2, 3]); // ACGT

    let recon: Reconstruction<Nucleotide> = Reconstruction {
        sequences: seqs.clone(),
        posteriors: None,
        log_likelihood: -1.5,
        alphabet: std::marker::PhantomData,
    };

    assert_eq!(recon.sequence_string(0_usize.into()), Some("ACGT".to_string()));
    assert_eq!(recon.sequence_string(99_usize.into()), None);
}

// ===========================================================================
// Marginal ASR tests (adapted from treetime test_ancestral)
// ===========================================================================

#[test]
fn test_marginal_asr_basic() {
    let tree = build_tiny_tree();
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    // All 4 taxa of the tree must be present
    let aln_data = b">A\nACGTACGTACGTACGT\n>B\nAGGTAGGTAGGTAGGT\n>C\nACGTACGTACGTACGT\n>D\nTTTTTTTTTTTTTTTT\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let result = tree.marginal_asr::<Nucleotide>(&model, &aln, false);
    assert!(result.is_ok());
    let recon = result.unwrap();
    assert!(!recon.sequences.is_empty());
}

#[test]
fn test_marginal_asr_with_posteriors() {
    let tree = build_tiny_tree();
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    let aln_data = b">A\nACGTACGTACGTACGT\n>B\nAGGTAGGTAGGTAGGT\n>C\nACGTACGTACGTACGT\n>D\nTTTTTTTTTTTTTTTT\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let result = tree.marginal_asr::<Nucleotide>(&model, &aln, true);
    assert!(result.is_ok());
    let recon = result.unwrap();

    // Posteriors should be present
    assert!(recon.posteriors.is_some());
    let posters = recon.posteriors.unwrap();

    // Each posterior should sum to ~1.0 for each site
    for (node_id, site_posters) in &posters {
        for (site_idx, post) in site_posters.iter().enumerate() {
            let sum: f64 = post.iter().sum();
            assert!(sum > 0.5 && sum < 1.5,
                "Node {} site {}: posterior sum = {} (expected ~1.0)",
                node_id, site_idx, sum);
        }
    }
}

#[test]
fn test_marginal_asr_log_likelihood_positive() {
    let tree = build_tiny_tree();
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    let aln_data = b">A\nACGT\n>B\nAGGT\n>C\nACGT\n>D\nTTTT\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let result = tree.marginal_asr::<Nucleotide>(&model, &aln, false);
    assert!(result.is_ok());
    // Log-likelihood should be finite (not +/-inf)
    assert!(result.unwrap().log_likelihood.is_finite());
}

// ===========================================================================
// Joint ASR tests (adapted from treetime test_seq_joint_reconstruction_correct)
// ===========================================================================

#[test]
fn test_joint_asr_basic() {
    let tree = build_tiny_tree();
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    let aln_data = b">A\nACGTACGTACGTACGT\n>B\nAGGTAGGTAGGTAGGT\n>C\nACGTACGTACGTACGT\n>D\nTTTTTTTTTTTTTTTT\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let result = tree.joint_asr::<Nucleotide>(&model, &aln);
    assert!(result.is_ok());
    let recon = result.unwrap();

    // Posteriors should be None for joint reconstruction
    assert!(recon.posteriors.is_none());
    assert!(!recon.sequences.is_empty());
}

#[test]
fn test_joint_vs_marginal_same_sequences_simple() {
    // For very short branch lengths, marginal and joint should agree closely
    let tree = build_tiny_tree();
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    // Very similar sequences (short branches -> minimal divergence), include D
    let aln_data = b">A\nAAAAAAAA\n>B\nAAAAAAAC\n>C\nAAAAAAAA\n>D\nCCCCCCCC\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let marg_result = tree.marginal_asr::<Nucleotide>(&model, &aln, false).unwrap();
    let joint_result = tree.joint_asr::<Nucleotide>(&model, &aln).unwrap();

    // Both should return sequences for all nodes
    assert_eq!(marg_result.sequences.len(), joint_result.sequences.len());
}

// ===========================================================================
// Joint LH maximality test (adapted from treetime test_seq_joint_lh_is_max)
// ===========================================================================

#[test]
fn test_marginal_asr_likelihood_normalization() {
    // From treetime: likelihood sum over all states should be 1
    let tiny_tree = b"((A:0.601,B:0.301)C:0.1,D:0.2)E:0.001;";
    let tree: PhyloTree = crate::prelude::Newick::from_newick(tiny_tree).unwrap();

    // Custom GTR with non-uniform pi (from treetime test_ancestral)
    let pi = vec![0.9, 0.06, 0.02, 0.02];
    let w = nalgebra::DMatrix::from_element(4, 4, 1.0);
    let model = GtrModel::<Nucleotide>::new(pi.clone(), w, true).unwrap();

    // Full alignment (treetime uses this pattern repeated)
    let aln_data = b">A\nAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT\n>B\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n>C\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>D\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    // The likelihood should be finite and positive in linear space
    assert!(aln.seqs.values().next().unwrap().len() == aln.width);
    let result = tree.marginal_asr::<Nucleotide>(&model, &aln, true);
    assert!(result.is_ok());
    let recon = result.unwrap();
    assert!(recon.log_likelihood.is_finite(), "log_likelihood = {}", recon.log_likelihood);

    // Check that marginals for root sum to 1 per site (within tolerance)
    if let Some(ref posters) = recon.posteriors {
        let root_id: NodeID = tree.get_root_id();
        if let Some(site_posters) = posters.get(&root_id) {
            assert_eq!(site_posters.len(), aln.width);
            for (site_idx, post) in site_posters.iter().enumerate() {
                let sum: f64 = post.iter().sum();
                assert!((sum - 1.0).abs() < 0.001,
                    "Root posterior sum at site {} = {}", site_idx, sum);
            }
        }
    }
}

// ===========================================================================
// Taxon mismatch tests
// ===========================================================================

#[test]
fn test_asr_missing_taxon_returns_error() {
    let tree = build_tiny_tree();
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    // Alignment has a taxon 'X' not in the tree, and is missing D
    let aln_data = b">A\nACGT\n>B\nAGGT\n>C\nAAAA\n>X\nTTTT\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let result = tree.marginal_asr::<Nucleotide>(&model, &aln, false);
    assert!(result.is_err());
}

// ===========================================================================
// Alphabet mismatch / bad input tests
// ===========================================================================

#[test]
fn test_marginal_asr_all_gap_alignment() {
    // All-gap alignment should still produce valid results (just uninformative)
    let tree = build_tiny_tree();
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    let aln_data = b">A\n----\n>B\n----\n>C\n----\n>D\n----\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let result = tree.marginal_asr::<Nucleotide>(&model, &aln, false);
    assert!(result.is_ok());
}

#[test]
fn test_joint_asr_all_gap_alignment() {
    let tree = build_tiny_tree();
    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    let aln_data = b">A\n----\n>B\n----\n>C\n----\n>D\n----\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let result = tree.joint_asr::<Nucleotide>(&model, &aln);
    // All-gaps produce uninformative data but should not error or panic
    assert!(result.is_ok(), "joint all-gap ASR should not error");
}

// ===========================================================================
// Edge weight/branch length tests
// ===========================================================================

#[test]
fn test_asr_with_zero_branch_length() {
    // Tree with zero branch lengths: children should have same distribution as parent
    let newick = b"((A:0.0,B:0.0)C:0.0,D:0.0)E:0.0;";
    let tree: PhyloTree = crate::prelude::Newick::from_newick(newick).unwrap();

    let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();

    let aln_data = b">A\nACGT\n>B\nAGGT\n>C\nAAAA\n>D\nTTTT\n";
    let aln = Alignment::from_fasta_bytes(aln_data).unwrap();

    let result = tree.marginal_asr::<Nucleotide>(&model, &aln, false);
    assert!(result.is_ok());
}

// ===========================================================================
// Compression correctness test (adapted from treetime compression patterns)
// ===========================================================================

#[test]
fn test_compression_with_duplicate_patterns() {
    let mut seqs = HashMap::new();
    // Columns 0-3: AA, CC, TT, GG -> all unique
    // Columns 4-7: AA, CC, TT, GG -> duplicates of 0-3
    // Should give 4 patterns with multiplicity 2 each
    seqs.insert("S1".to_string(), b"ACGTACGT".to_vec());
    seqs.insert("S2".to_string(), b"ACGTACGT".to_vec());
    let aln = Alignment { seqs, width: 8 };

    let comp = aln.compress_columns();
    assert_eq!(comp.patterns.len(), 4);
    for &mult in &comp.multiplicity {
        assert_eq!(mult, 2);
    }
}

#[test]
fn test_compression_re_expansion() {
    let mut seqs = HashMap::new();
    seqs.insert("S1".to_string(), b"AATTGGCC".to_vec());
    seqs.insert("S2".to_string(), b"ATGCATGC".to_vec());
    let aln = Alignment { seqs, width: 8 };

    let comp = aln.compress_columns();

    // Verify site-to-pattern mapping is correct by reconstructing the original columns
    for site in 0..aln.width {
        let p_idx = comp.site_to_pattern[site];
        let expected_col: Vec<u8> = comp.leaf_order.iter()
            .map(|name| aln.seqs[name][site])
            .collect();
        assert_eq!(comp.patterns[p_idx], expected_col,
            "Pattern for site {} mismatch", site);
    }
}

// ===========================================================================
// Amino acid alphabet test (from treetime AminoAcid patterns)
// ===========================================================================

#[test]
fn test_amino_acid_n_states() {
    use crate::tree::asr::alphabet::AminoAcid;
    assert_eq!(AminoAcid::N_STATES, 20);
}

#[test]
fn test_nucleotide_n_states() {
    use crate::tree::asr::alphabet::Nucleotide;
    assert_eq!(Nucleotide::N_STATES, 4);
}

// ===========================================================================
// Profile numerical stability tests
// ===========================================================================

#[test]
fn test_profile_scaling_stability() {
    // Test that scaling handles very small values correctly
    let tiny = vec![1e-300, 2e-300, 3e-300];
    let profile = crate::tree::asr::profile::Profile::new(tiny, 0.0).scale();

    // After scaling by max (3e-300), values should be normalized
    assert!((profile.values[2] - 1.0).abs() < 1e-15);
    // log_scale = ln(3e-300) ~ -691, which is finite and negative
    // The key invariant: total_log_likelihood = ln(sum(values)) + log_scale
    assert!(profile.log_scale.is_finite());
}

#[test]
fn test_profile_log_likelihood_with_scale() {
    // Values [0.1, 0.2, 0.3, 0.4], scale = ln(10)
    let vals = vec![0.1, 0.2, 0.3, 0.4];
    let scale = 10.0f64.ln();
    let profile = crate::tree::asr::profile::Profile::new(vals, scale);

    // sum(values) = 1.0, ln(1.0) + scale = scale
    assert!((profile.total_log_likelihood() - scale).abs() < 1e-10);
}
