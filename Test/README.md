# Test dataset

The following genomes were used to create the DESwoMAN dataset:
- Queries:  _Drosophila melanogaster_: FI and ZI genomes and transcriptomes ([downloaded from Grandchamp et al. (2023)](https://zenodo.org/records/7322757))  
- Targets: 
    - _Drosophila simulans_(GCF_016746395.2)
    - _Drosophila suzukii_(GCF_037355615.1)
    - _Drosophila ananassae_(GCF_017639315.1)

DESwoMAN was run using the following parameters:
- Intergenic ORFs only
- TPM >= 0.5 (+ Exclude undirectional transcripts ".")
- Protein search (--more-sensitive option) against _Drosophila_ and Diptera from Ensembl Metazoa (database not uploaded due to size, available upon request)
- Synteny window = 4 (no reciprocal blast)

Additional data downloaded:
- [Flybase TE database](https://flybase.org/downloads/bulkdata) (accessed July, 2025)
- ncRNA database: [Drosophila simulans ncRNA](https://metazoa.ensembl.org/Drosophila_simulans_gca016746395v2rs/Info/Index) (Assembly: GCF_016746395.2)
