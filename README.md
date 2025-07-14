# IPFE Simulation

This directory contains a Python-based simulation of **Inner Product Functional Encryption (IPFE)**.  
It is designed for learning and experimentation purposes, implementing core concepts such as encryption, decryption, and discrete logarithm recovery using multiple strategies (Baby-Step Giant-Step, Pollard's Kangaroo, etc.).

> **Note:** This is a simulation — not a production-ready implementation. All parameters used (including the modulus, generator, and vector bounds) are small and for demonstration only.

### 🔧 Parameters
- `p` is a small prime modulus (e.g. 1009, 1019, etc.)
- `g` is a generator for the multiplicative group ℤₚ*
- `ptx` is the plaintext vector (e.g. `[3, 2, 1]`)
- `y` is the function vector for which the decryption yields ⟨x, y⟩
- `r` is random for encryption randomness
- All inner products are kept **below `p`** to ensure correctness of decryption

### 🧪 Features
- Key generation (`msk`, `mpk`)
- Functional encryption of vectors
- Recovery of ⟨x, y⟩ from ciphertext using known secret key
- Discrete log recovery via fallback strategy (brute force → BSGS → Pollard)

### 📚 Further Reading
To understand the theoretical basis of IPFE, refer to:

**Abdalla, M., et al. (2015).**  
*Simple Functional Encryption Schemes for Inner Products*  
[https://eprint.iacr.org/2015/017](https://eprint.iacr.org/2015/017)

---

This simulation is meant for educational use and experimentation only.
