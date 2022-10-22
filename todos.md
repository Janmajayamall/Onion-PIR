To understand

1. Understand how CBD function here works - https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-util/src/lib.rs#L37-L38

2. How does this work - https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/mod.rs#L426 ?

3. Why de we add `self.params.plaintext_modulus_u64` to `a` (before reducing resulting value by `mod q[0]`) when converting descaled plaintext from `mod q[0]` to `mod t`?

4. Understand implementation of U256, and `fn scale` for `Scalar`.

5. Understand how moduli are generated for moduli bit sizes.

6. Understand sample_vec_cbd

7. Understand substitution in NTT domain.

8. Why does key switching (and substitution) require ciphertext modulues residue moduli length be greater than 1 ?
   Is it because error growth is too much for ciphertext space with 64 bit to handle ?
   Ref - https://github.com/tlepoint/fhe.rs/blob/f7cddb358f2ce28483944f99e223c07ae41b0c1c/crates/fhe/src/bfv/keys/key_switching_key.rs#L58

9.

Work on

1. Implement substitutions and unpack algorithm of onion PIR
2. Understand (4), (7), (1,6)
