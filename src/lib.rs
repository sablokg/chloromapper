//! `chloroplast_ir`: given the assembled FASTA of a chloroplast genome,
//! locate the inverted-repeat (IRa/IRb) regions and emit a single,
//! canonically-ordered, end-to-end genome sequence:
//!
//! ```text
//! LSC -> IRb -> SSC -> IRa -> (back to LSC, circular)
//! ```
//!
//! Two input shapes are supported:
//!
//! * **One already-complete circular contig** (the common case for modern
//!   organelle assemblers such as GetOrganelle/NOVOPlasty/Unicycler): the
//!   genome is simply *rotated* to start at the canonical LSC boundary — no
//!   base is invented and no `N` is needed, because a correctly assembled
//!   single contig is already self-consistent on one strand.
//!
//! * **Several separate contigs** (e.g. the assembly graph split the
//!   genome at repeat boundaries): contigs are classified into LSC / IR /
//!   SSC roles, stitched into canonical order, and `N` is inserted at each
//!   inter-contig junction (configurable, default 0) since the exact
//!   junction sequence is not recoverable from the contigs alone. If only
//!   one IR copy was assembled (common — the graph often collapses the
//!   repeat into a single contig), the second copy is synthesized as its
//!   exact reverse complement.

pub mod arrange;
pub mod classify;
pub mod fasta;
pub mod ir_finder;
pub mod seq;

/*
Gaurav Sablok
gsablok@proton.me
*/

pub use arrange::{arrange_from_contigs, arrange_single_contig, ArrangedGenome, RegionMap};
pub use classify::{classify, Classified};
pub use fasta::{FastaError, FastaRecord};
pub use ir_finder::{find_inverted_repeats, IrFinderConfig, IrPair};
