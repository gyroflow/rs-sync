// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright © 2021-2022 Vladimir Pinchuk (https://github.com/VladimirP1)
// Ported to Rust by Adrian <adrian.eddy at gmail>

#![allow(unused)]

use nalgebra::*;

#[derive(Default)]
pub struct Spline {
    m_y: Vec<f64>,
    m_b: Vec<f64>,
    m_c: Vec<f64>,
    m_d: Vec<f64>
}

impl Spline {
    pub fn new(y: &[f64]) -> Self {
        let mut ret = Self::default();
        ret.set_points(y);
        ret
    }
    fn set_points(&mut self, y: &[f64]) {
        let n = y.len();
        if n == 0 { return; }

        let mut a = MatrixXx3::<f64>::from_element(n, 0.0);
        self.m_c.resize(n, 0.0);
        for i in 1..n - 1 {
            a[(i, 0)] = 1.0 / 3.0;
            a[(i, 1)] = 2.0 / 3.0 * 2.0;
            a[(i, 2)] = 1.0 / 3.0;
            self.m_c[i] = y[i + 1] - 2.0 * y[i] + y[i - 1];
        }

        a[(0, 1)] = 2.0;
        a[(0, 2)] = 0.0;
        self.m_c[0] = 0.0;

        a[(n - 1, 1)] = 2.0;
        a[(n - 1, 0)] = 0.0;
        self.m_c[n - 1] = 0.0;

        for i in 0..n - 2 {
            let k = 1.0 / a[(i, 1)] * a[(i + 1, 0)];
            a[(i + 1, 0)] -= a[(i, 1)] * k;
            self.m_c[i + 1] -= self.m_c[i] * k;
        }

        for i in (1..n - 1).rev() {
            let k = 1.0 / a[(i, 1)] * a[(i - 1, 2)];
            a[(i - 1, 1)] -= a[(i, 0)] * k;
            self.m_c[i - 1] -= self.m_c[i] * k;
        }

        for (i, x) in self.m_c.iter_mut().enumerate() {
            *x /= a[(i, 1)];
        }

        self.m_d.resize(n, 0.0);
        self.m_b.resize(n, 0.0);
        for i in 0..n - 1 {
            self.m_d[i] = 1.0 / 3.0 * (self.m_c[i + 1] - self.m_c[i]);
            self.m_b[i] = (y[i + 1] - y[i]) - 1.0 / 3.0 * (2.0 * self.m_c[i] + self.m_c[i + 1]);
        }

        self.m_d[n - 1] = 0.0;
        self.m_b[n - 1] = 3.0 * self.m_d[n - 2] + 2.0 * self.m_c[n - 2] + self.m_b[n - 2];
        self.m_y = y.to_vec();
    }

    pub fn eval(&self, x: f64) -> Option<f64> {
        let idx = x.floor().min(self.m_b.len() as f64).max(0.0) as usize;
        let n = self.m_b.len();
        let h = x - idx as f64;
        if x < idx as f64 { return Some((self.m_c.first()? * h + self.m_b.first()?) * h + self.m_y.first()?); }
        if x > n as f64 - 1.0 { return Some((self.m_c.get(n - 1)? * h + self.m_b.get(n - 1)?) * h + self.m_y.get(n - 1)?); }
        Some(((self.m_d.get(idx)? * h + self.m_c.get(idx)?) * h + self.m_b.get(idx)?) * h + self.m_y.get(idx)?)
    }
    pub fn deriv(&self, x: f64) -> Option<f64> {
        let idx = x.floor().min(self.m_b.len() as f64).max(0.0) as usize;
        let n = self.m_b.len();
        let h = x - idx as f64;
        if x < 0.0 { return Some(2.0 * self.m_c.first()? * h + self.m_b.first()?); }
        if x > n as f64 - 1.0 { return Some(2.0 * self.m_c.get(n - 1)? * h + self.m_b.get(n - 1)?); }
        Some((3.0 * self.m_d.get(idx)? * h + 2.0 * self.m_c.get(idx)?) * h + self.m_b.get(idx)?)
    }
}
