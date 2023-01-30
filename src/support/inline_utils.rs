// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright © 2021-2022 Vladimir Pinchuk (https://github.com/VladimirP1)
// Ported to Rust by Adrian <adrian.eddy at gmail>

#![allow(unused)]

use nalgebra::*;
use rand::prelude::*;

pub fn safe_normalize(m: Vector3<f64>) -> Vector3<f64> {
    let norm = m.norm();
    if norm < 1e-12 {
        return m;
    }
    m / norm
}

pub fn mtrand(min: usize, max: usize) -> usize {
    let mut rng = thread_rng();
    let distr = rand::distributions::Uniform::new_inclusive(min, max);
    rng.sample(distr)
}

pub fn sqr_jac(x: &DVector<f64>) -> (DVector<f64>, SquareMatrix<f64, Dyn, VecStorage<f64, Dyn, Dyn>>) {
    (x.map(|v| v * v), SquareMatrix::from_diagonal(&(2.0 * x)))
}

pub fn sqrt_jac(mut x: DVector<f64>) -> (DVector<f64>, SquareMatrix<f64, Dyn, VecStorage<f64, Dyn, Dyn>>) {
    x.apply(|v| *v = v.sqrt());
    let diag = x.map(|v| 1.0 / (2.0 * v));
    (x, SquareMatrix::from_diagonal(&diag))
}

pub fn log1p_jac(x: DVector<f64>) -> (DVector<f64>, SquareMatrix<f64, Dyn, VecStorage<f64, Dyn, Dyn>>) {
    (x.map(libm::log1p), SquareMatrix::from_diagonal(&x.map(|v| 1.0 / (1.0 + v))))
}

pub fn sum_jac(x: &DMatrix<f64>) -> (Matrix1xX<f64>, Matrix1xX<f64>) {
    (x.row_sum(), Matrix1xX::from_element(x.nrows(), 1.0))
}

pub fn div_jac(x: &DMatrix<f64>, y: f64) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
    let dx = DMatrix::identity(x.nrows(), x.nrows());
    (x / y, dx / y, -x / (y * y))
}

pub fn mul_const_jac(x: &DMatrix<f64>, y: f64) -> (DMatrix<f64>, DMatrix<f64>) {
    let dx = DMatrix::identity(x.nrows(), x.nrows());
    (x * y, dx * y)
}

pub fn clamp_k(k : f64) -> f64 {
    k.clamp(1e1, 1e3)
}