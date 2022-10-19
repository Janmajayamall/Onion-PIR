To understand

1. Understand how CBD function here works - https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-util/src/lib.rs#L37-L38

2. How does this work - https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/mod.rs#L426 ?

3. Why de we add `self.params.plaintext_modulus_u64` to `a` (before reducing resulting value by `mod q[0]`) when converting descaled plaintext from `mod q[0]` to `mod t`?

4. Understand implementation of U256, and `fn scale` for `Scalar`.

5. Understand how moduli are generated for moduli bit sizes.

Work on

1. Plaintext (and what are levels)
