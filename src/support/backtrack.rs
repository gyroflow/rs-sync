// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright © 2021-2022 Vladimir Pinchuk (https://github.com/VladimirP1)
// Ported to Rust by Adrian <adrian.eddy at gmail>

struct BacktrackHyper {
    sufficent_decrease: f64,
    decay: f64,
    initial_step: f64,
    max_iterations: usize
}
impl Default for BacktrackHyper {
    fn default() -> Self {
        Self {
            sufficent_decrease: 0.7,
            decay: 0.1,
            initial_step: 1.0,
            max_iterations: 20
        }
    }
}
#[derive(Default)]
pub struct Backtrack<'a> {
    hyper: BacktrackHyper,
    f_and_grad: Option<Box<dyn Fn(f64) -> (f64, f64) + 'a>>,
    f_only: Option<Box<dyn Fn(f64) -> f64 + 'a>>,
}

impl<'a> Backtrack<'a> {
    pub fn set_hyper(&mut self, sufficent_decrease: f64, decay: f64, initial_step: f64, max_iterations: usize) {
        self.hyper.sufficent_decrease = sufficent_decrease;
        self.hyper.decay = decay;
        self.hyper.initial_step = initial_step;
        self.hyper.max_iterations = max_iterations;
    }

    pub fn set_objective_f_only<F>(&mut self, f_only: F) where F: Fn(f64) -> f64 + 'a { self.f_only = Some(Box::new(f_only)); }

    pub fn set_objective<F>(&mut self, f_and_grad: F) where F: Fn(f64) -> (f64, f64) + Clone + 'a {
        let f_and_grad2 = f_and_grad.clone();
        self.f_and_grad = Some(Box::new(f_and_grad));
        if self.f_only.is_none() {
            self.f_only = Some(Box::new(move |x: f64| { f_and_grad2(x).0 }));
        }
    }

    pub fn step(&self, x0: f64) -> f64 {
        let (v, p) = self.f_and_grad.as_ref().unwrap()(x0);
        let m = p * p;
        let mut t = self.hyper.initial_step;
        for _ in 0..self.hyper.max_iterations {
            let v1 = self.f_only.as_ref().unwrap()(x0 - t * p);
            if v - v1 >= t * self.hyper.sufficent_decrease * m { break; }
            t *= self.hyper.decay;
        }
        -t * p
    }
}
