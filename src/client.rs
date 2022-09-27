use std::sync::Arc;

use crate::{
    bfv::{BfvCipherText, BfvParameters, BfvPlaintext, BfvPublicKey},
    poly::Poly,
};

struct QueryParams {
    first_dim: u64,
    normal_dim: u64,
}

fn build_query(
    query_params: QueryParams,
    mut db_dim: u64,
    mut query_index: u64,
    bfv_pk: &BfvPublicKey,
) -> BfvCipherText {
    let mut no_of_dims = 0;
    let mut first = true;
    while db_dim != 0 {
        if first {
            db_dim /= query_params.first_dim;
            first = false;
        } else {
            db_dim /= query_params.normal_dim;
        }
        no_of_dims += 1;
    }

    let mut dims = Vec::new();
    (1..no_of_dims + 1).into_iter().for_each(|dim| {
        if dim == 1 {
            dims.push(query_index % query_params.first_dim);
            query_index /= query_params.first_dim;
        } else {
            dims.push(query_index % query_params.normal_dim);
            query_index /= query_params.normal_dim
        }
    });

    // Imagine the db vector as a hypercube of `d` dimensions:
    // [d0, d1, d2, d3, d4...d{d-1}].
    // To reach a specific query index you iteratively
    // divide current db dimension with next hypercube dimension.
    // At each iteration (i.e. for current dimension) your dimension specific query
    // index is equal to last query index `mod` current db dimension.
    // You next query index is equal to last query index / current db dimension.
    //
    // Once you have obtained vector containing dimension specific query
    // indexes. You need to expand them into a vec to 0s and 1s for bfv ciphertext
    // encoding. For example, consider following as dimension specific query indexes
    // [8, 2, 3].
    // You expand them to
    // [
    //     [
    //         // first_dim len
    //         0,0,0,0,0,0,0,1,0,0,...,0
    //     ],
    //     [
    //         // normal_dim len
    //         0,0,1,0
    //     ],
    //     [
    //         // normal_dim len
    //         0,0,0,1
    //     ]
    // ]
    // and encode them into Bfv ciphertext
    let mut expanded_dims = Vec::new();
    dims.iter().enumerate().for_each(|(index, value)| {
        let mut expanded_dim = Vec::<u64>::new();
        (0..(if index == 0 {
            query_params.first_dim
        } else {
            query_params.normal_dim
        }))
            .into_iter()
            .for_each(|i| {
                if i == *value {
                    expanded_dim.push(1);
                } else {
                    expanded_dim.push(0);
                }
            });
        expanded_dims.push(expanded_dim);
    });

    // encode dimension specific query vector bits into bfv ciphertext
    let pt = BfvPlaintext::new(
        &bfv_pk.params,
        expanded_dims.iter().flatten().copied().collect(),
    );

    bfv_pk.encrypt(&pt)
}
