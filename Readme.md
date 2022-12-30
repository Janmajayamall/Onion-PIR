# ONION-PIR

This repository implements prototype of [Onion-PIR]() from scratch. This means along with algorithms used in OnionPIR, you will find implementation of primitives that PIR schemes use in general.

To test checkout `evaluations.rs`, more specifically `evaluations::tests::query`.

## Oblivious expansion

```rb

Let's consider the following bit vector
[0, 1, 0, 1, 1, 0...]

 We encode the bit vector into a `BfvPlaintext` as following
 `0•X^0 + 1•X^1 + 0•X^2 + 1•X^3 + ... + Nth_bit•X^{N-1}
 and then encrypt it under secret key `SK`. Let's call the
 encrypted ciphertext `C`.

 We need to recover individual bits encrypted under `SK`
 from the ciphertext `C`. We achieve this using the following
 way:

 Assume that our poly operations are modulo cyclotomic
 polynomial `X^16 + 1`. (i.e. N = 16)

 Note: that `X^N = -1` and `(X^i)^N+1 = X^(Ni+i) = X^Ni * X^i`
 Therefore, `(X^i)^N+1 = -1^i * X^i`

 We will use the statement above to separate encryption of odd
 terms and even terms in `C` into two different encryptions.

 `C_even = C + Subs(C, N+1)` --> Encryption of only even terms
 `C_odd = C - Subs(C, N+1)` --> Encryption of only odd terms

 Note: that `Subs(C, N+1)` transforms `Enc(Σ b_i • X^i)` to
 `Enc(Σ b_i • (X^i)^(N+1))`. Also, `Σ b_i • (X^i)^(N+1)` is equal
 to `Σ b_i • (X^i)` except that all its odd position are -ve.
 Therefore,
 `Enc(Σ b_i • X^i)` + `Enc(Σ b_i • (X^i)^(N+1))` = Ct_even
 AND
 `Enc(Σ b_i • X^i)` - `Enc(Σ b_i • (X^i)^(N+1))` = Ct_odd

 Also notice that in case of `Ct_odd` the powers are odd, which
 wouldn't isn't useful in subsequent steps. Thus we change Ct_odd
 `Ct_odd = Ct_odd * X^-1`

 Therefore now we have
 `Ct_even = C + Subs(C, N+1) = Enc(Σ 2b_2i • X^2i)`
 `Ct_odd = (C - Subs(C, N+1)) • X^-1 = Enc(Σ 2b_2i+1 • X^2i)`

 Now to obtain `Enc(b_i)`s of bits we perform the same operation
 recursively. After `log(N)` branches we will have `Enc(N * b_i)`
 as N leaf nodes.
 `k = n/2^(i) + 1`
                                    [x0, x1, x2, x3, x4, x5, x6, x7] depth of tree: logN
                                      / i = 0                                       \
                                     / C + Subs(•, k)                                \ C - Subs(•, k)
                                  2^(i+1)*[x0, x2, x4, x6]                     2^(i+1)*[x1, x3, x5, x7] * x^-(2^i) => 2^i*[x0, x2, x4, x6]
                                    / i = 1              \
                                   / C + Subs(•, k)       \ C - Subs(•, k)
                             2^(i+1)*[x0, x4]             2^(i+1)*[x2, x6] * x^-(2^i)
                             / i = 2        \
                            / C + Subs(•, k) \ C - Subs(•, k)
                     2^(i+1)*[x0]         2^(i+1)*[x4] * x^-(2^i)

 Ref - Algorithm 3 & 4 of https://eprint.iacr.org/2019/736.pdf
```

## How does query work?

```rb
/// Imagine the db vector as a hypercube of `d` dimensions:
/// [d0, d1, d2, d3, d4...d{d-1}].
/// To reach a specific query index you iteratively
/// divide current db dimension with next hypercube dimension.
/// At each iteration (i.e. for current dimension) your dimension specific query
/// index is equal to last query index `mod` current db dimension.
/// You next query index is equal to last query index / current db dimension.
///
/// Once you have obtained vector containing dimension specific query
/// indexes. You need to expand them into a vec to 0s and 1s for bfv ciphertext
/// encoding. For example, consider following as dimension specific query indexes
/// [8, 2, 3].
/// You expand them to
/// [
///     [
///         // first_dim len
///         0,0,0,0,0,0,0,1,0,0,...,0
///     ],
///     [
///         // normal_dim len
///         0,0,1,0
///     ],
///     [
///         // normal_dim len
///         0,0,0,1
///     ]
/// ]
/// and encode them into Bfv ciphertext. This means simply flat
/// map all values as coefficients of the Bfv plaintext
/// polynomial
```

## Credits

1. [fhe.rs](https://github.com/tlepoint/fhe.rs): Code is heavily inspired from fhe.rs.
2. [lattigo](https://github.com/tuneinsight/lattigo): Several implementation details were inspired from lattigo.
