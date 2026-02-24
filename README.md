# QBastion (Quantum Bastion)

Quantum Node Protection system using STARK-based aggregation of quantum-resistant signatures.

This project adapts the s2morrow architecture to create a quantum-resistant node system. The goal is to detect and block quantum computers attempting to decode cryptocurrency wallets, utilizing zkVM for batch verification of multiple PQ signatures. 

Implementation details:
- ZKVM: Cairo
- STARK prover: Stwo
- Signature schemes: NIST finalist — Falcon (lattice based), alternative Sphincs+ (hash based)

## Roadmap

- [x] Falcon512 verification
- [x] Sphincs+ 128s with SHA2 (simple mode) verification
- [ ] Stwo proving benchmarks

Follow-up:
- [x] Sphincs+ 128s with Blake2s and 4-byte aligned address encoding
- [ ] Falcon512 with probabilistic polynomial multiplication checking

## References

- [BIP360](https://bip360.org/) Pay to Quantum Resistant Hash
- [PQ Signatures and Scaling Bitcoin with STARKs](https://delvingbitcoin.org/t/post-quantum-signatures-and-scaling-bitcoin-with-starks/1584)
- Related [thread](https://groups.google.com/g/bitcoindev/c/wKizvPUfO7w/m/hG9cwpOABQAJ) in bitcoindev mailing list
- Hash based PQ signatures [part 1](https://research.dorahacks.io/2022/10/26/hash-based-post-quantum-signatures-1/) and [part 2](https://research.dorahacks.io/2022/12/16/hash-based-post-quantum-signatures-2/)