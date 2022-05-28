// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright © 2021-2022 Vladimir Pinchuk (https://github.com/VladimirP1)
// Ported to Rust by Adrian <adrian.eddy at gmail>

#![allow(unused)]

use super::minispline::*;
use super::quat::*;
use nalgebra::*;

#[derive(Default)]
pub struct NdSpline {
    splines: Vec<Spline>
}

impl NdSpline {
    pub fn make(m: &Matrix4xX<f64>) -> Self {
        let mut ret = Self::default();
        for row in m.row_iter() {
            ret.splines.push(Spline::new(row.transpose().as_slice()))
        }
        dbg!(ret.splines.len());
        ret
    }
    pub fn eval(&self, t: f64) -> Vector4<f64> {
        let mut ret = Vector4::from_element(0.0);
        for (i, sp) in self.splines.iter().enumerate() {
            ret[i] = sp.call(t);
        }
        ret
    }
    pub fn deriv(&self, t: f64) -> Vector4<f64> {
        let mut ret = Vector4::from_element(0.0);
        for (i, sp) in self.splines.iter().enumerate() {
            ret[i] = sp.deriv(t);
        }
        ret
    }
    pub fn rderiv_numeric(&self, t: f64) -> Vector4<f64> {
        let i_l = Vector4::from(self.eval(t).normalize());
        let i_r = Vector4::from(self.eval(t + 1e-7).normalize());
        let mut ret = quat_prod(&quat_conj(i_l), &i_r) / 1e-7;
        ret[0] = 0.0;
        ret
    }
    pub fn rderiv(&self, t: f64) -> Vector4<f64> {
        let value = self.eval(t);
        let norm = value.norm();
        (quat_prod(&quat_conj(value), &self.deriv(t))) / (norm * norm) * 2.0
    }
}
