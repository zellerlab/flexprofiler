# PathSeq-T2T resources

`non_human_decoys.bed` lists the viral decoy contigs of the GRCh38 reference of the NCI Genomic Data Commons
(GRCh38.d1.vd1). Reads aligned to them are not host reads, and are given back to the filtering in the PathSeq-T2T mode
(`--pathseq_t2t`, see `docs/pathseq_t2t.md`).

It is copied unmodified from [PathSeq-T2T](https://github.com/abdohlman/pathseq-t2t) (commit
`51d74430b6a4b34073f2d878612ac24b6a1d1e80`), which is distributed under the MIT license in `LICENSE.pathseq-t2t`.

> Dohlman A.B., Mjelle R., Wood H.M., et al. Biodiversity and biogeography of the multi-kingdom cancer microbiome.
> Cell 189, 3817-3837 (2026). doi: [10.1016/j.cell.2026.04.015](https://doi.org/10.1016/j.cell.2026.04.015)
