// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright © 2021-2022 Vladimir Pinchuk (https://github.com/VladimirP1)
// Ported to Rust by Adrian <adrian.eddy at gmail>

mod support {
    pub(crate) mod backtrack;
    pub(crate) mod minispline;
    pub(crate) mod ndspline;
    pub(crate) mod quat;
    pub(crate) mod inline_utils;
}
use support::{ backtrack::*, ndspline::*, quat::*, inline_utils::* };

use nalgebra::*;
use superslice::*;
use argmin::{ core::{ Gradient, CostFunction, Executor, Error }, solver::linesearch::*, solver::linesearch::condition::ArmijoCondition, solver::quasinewton::LBFGS };
use std::collections::BTreeMap;
use std::cell::RefCell;
use rayon::iter::*;

pub struct FrameData {
    ts_a: DVector<f64>,
    ts_b: DVector<f64>,
    rays_a: Vec<Vector3<f64>>,
    rays_b: Vec<Vector3<f64>>,
}

#[derive(Default)]
pub struct OptData {
    quats_start: f64,
    sample_rate: f64,
    quats: NdSpline,

    frame_data: BTreeMap<i64, FrameData>
}

const NUMERIC_DIFF_STEP: f64 = 1e-6;

pub struct FrameState<'a> {
    timestamp_us: i64,
    problem: &'a OptData,

    motion_vec: DVector<f64>,
    var_k: f64,
}

impl<'a> FrameState<'a> {
    pub fn new(timestamp_us: i64, problem: &'a OptData) -> Self {
        Self {
            timestamp_us,
            problem,
            motion_vec: DVector::from_element(3, 0.0),
            var_k: 1e3
        }
    }

    pub fn loss(&self, gyro_delay: f64, motion_estimate: &DVector<f64>, jac_gyro_delay: &mut f64, jac_motion_estimate: &mut DVector<f64>) -> f64 {
        let p = opt_compute_problem(self.timestamp_us, gyro_delay, &self.problem);

        let loss_l = self.loss_single(gyro_delay - NUMERIC_DIFF_STEP, motion_estimate);
        let loss_r = self.loss_single(gyro_delay + NUMERIC_DIFF_STEP, motion_estimate);

        let (v1, j1) = (&p * motion_estimate, p);
        let (v2, j2) = sqr_jac(&v1);

        let (v3, j3) = sqr_jac(motion_estimate);
        let (v4, j4) = sum_jac(&convert(v3));
        let (v5, j5, _) = div_jac(&convert(v4), self.var_k * self.var_k);

        let (v6, j6a, j6b) = div_jac(&convert(v2), v5[0]);
        let (v7, j7) = log1p_jac(convert(v6));
        let (v8, j8) = sum_jac(&convert(v7));

        *jac_gyro_delay = (loss_r - loss_l) / 2.0 / NUMERIC_DIFF_STEP;

        *jac_motion_estimate = convert((j8 * j7 * (j6a * j2 * j1 + j6b * j5 * j4 * j3)).transpose());

        v8[0]
    }

    pub fn loss_single(&self, gyro_delay: f64, motion_estimate: &DVector<f64>) -> f64 {
        let p = opt_compute_problem(self.timestamp_us, gyro_delay, &self.problem);
        let r = (p * motion_estimate) * (self.var_k / motion_estimate.norm());
        let rho = r.map(|v| libm::log1p(v * v));
        rho.sum()
    }

    pub fn guess_motion(&self, gyro_delay: f64) -> DVector<f64> {
        let p = opt_compute_problem(self.timestamp_us, gyro_delay, &self.problem);
        opt_guess_translational_motion(&p, 200)
    }

    pub fn guess_k(&self, gyro_delay: f64) -> f64 {
        let p = opt_compute_problem(self.timestamp_us, gyro_delay, &self.problem);
        assert!(p.nrows() > 0);
        assert!(self.motion_vec.nrows() == 3);
        clamp_k(1.0 / (p * &self.motion_vec).norm() * 1e2)
    }

}

#[derive(Default)]
pub struct SyncProblem<'a> {
    problem: OptData,
    progress_cb: Option<Box<dyn Fn(f64) -> bool + Sync + 'a>>
}

impl<'a> SyncProblem<'a> {
    pub fn new() -> Self { Self::default() }
    pub fn on_progress<F: Fn(f64) -> bool + Sync + 'a>(&mut self, cb: F) {
        self.progress_cb = Some(Box::new(cb));
    }
    pub fn set_gyro_quaternions_fixed(&mut self, data: &[(f64, f64, f64, f64)], sample_rate: f64, first_timestamp: f64) {
        let flat_data = data.into_iter().flat_map(|x| [x.0, x.1, x.2, x.3]).collect::<Vec<f64>>();
        self.problem.sample_rate = sample_rate;
        self.problem.quats_start = first_timestamp;
        self.problem.quats = NdSpline::make(&Matrix4xX::from_column_slice(&flat_data));
    }

    pub fn set_gyro_quaternions(&mut self, timestamps_us: &[i64], quats: &[(f64, f64, f64, f64)]) {
        if timestamps_us.is_empty() {
            log::error!("Empty timestamps! quats.len: {}", quats.len());
            return;
        }

        const UHZ_IN_HZ: i64 = 1000000;
        const US_IN_SEC: i64 = 1000000;
        let count = timestamps_us.len();

        let actual_sr_uhz = UHZ_IN_HZ * US_IN_SEC * count as i64 / (timestamps_us[count - 1] - timestamps_us[0]);
        let rounded_sr_hz = ((actual_sr_uhz as f64 / 50.0 / UHZ_IN_HZ as f64).round() * 50.0) as i64;  // round to nearest 50hz

        if rounded_sr_hz <= 0 {
            log::error!("Invalid sample rate, count: {}, ts diff: {}", count, (timestamps_us[count - 1] - timestamps_us[0]));
            return;
        }

        let mut new_timestamps_vec = Vec::new();
        let mut sample = timestamps_us[0] * rounded_sr_hz / US_IN_SEC;
        while US_IN_SEC * sample / rounded_sr_hz < timestamps_us[count - 1] {
            new_timestamps_vec.push(US_IN_SEC * sample / rounded_sr_hz);
            sample += 1;
        }

        for i in 1..count {
            if timestamps_us[i - 1] > timestamps_us[i] {
                log::error!("timestamps out of order at pos {i} ({} > {})", timestamps_us[i - 1], timestamps_us[i]);
            }
        }

        let mut new_quats = Matrix4xX::from_element(new_timestamps_vec.len(), 0.0);
        for (i, ts) in new_timestamps_vec.iter().enumerate() {
            let idx = timestamps_us.lower_bound(ts);
            if idx > 0 {
                let t = 1.0 * (ts - timestamps_us[idx - 1]) as f64 / (timestamps_us[idx] - timestamps_us[idx - 1]) as f64;
                let a = Vector4::new(quats[idx - 1].0, quats[idx - 1].1, quats[idx - 1].2, quats[idx - 1].3);
                let b = Vector4::new(quats[idx].0, quats[idx].1, quats[idx].2, quats[idx].3);

                new_quats.set_column(i, &quat_slerp(&a, b, t));
            } else {
                new_quats.set_column(i, &Vector4::new(quats[idx].0, quats[idx].1, quats[idx].2, quats[idx].3));
            }
            // panic_to_file("set-gyro-quaternions: non-finite sample after interpolation", !new_quats.col(i).is_finite());
        }
        if new_timestamps_vec.is_empty() {
            log::error!("Invalid new timestamps: first: {}, last: {}, len: {}", timestamps_us[0], timestamps_us[count - 1], count);
            return;
        }
        self.problem.sample_rate = 1.0 * rounded_sr_hz as f64;
        self.problem.quats_start = 1.0 * new_timestamps_vec[0] as f64 / US_IN_SEC as f64;
        // panic_to_file("set-gyro-quaternions: non-finite sample rate. wtf?", !std::isfinite(problem.sample_rate));
        // panic_to_file("set-gyro-quaternions: non-finite first timestamp. wtf?", !std::isfinite(problem.quats_start));
        self.problem.quats = NdSpline::make(&new_quats);
    }

    pub fn set_track_result(&mut self, timestamp_us: i64, ts_a: &[f64], ts_b: &[f64], rays_a: &[(f64, f64, f64)], rays_b: &[(f64, f64, f64)]) {
        assert!(ts_a.len() == ts_b.len());
        assert!(rays_a.len() == rays_b.len());

        self.problem.frame_data.insert(timestamp_us, FrameData {
            rays_a:  rays_a.iter().map(|&(x, y, z)| Vector3::new(x, y, z)).collect(),
            rays_b:  rays_b.iter().map(|&(x, y, z)| Vector3::new(x, y, z)).collect(),
            ts_a:    DVector::from_row_slice(ts_a),
            ts_b:    DVector::from_row_slice(ts_b)
        });
        // panic_to_file("set-track-result: non-finite numbers in rays_a", !flow.rays_a.is_finite());
        // panic_to_file("set-track-result: non-finite numbers in rays_b", !flow.rays_b.is_finite());
        // panic_to_file("set-track-result: non-finite numbers in ts_a", !flow.ts_a.is_finite());
        // panic_to_file("set-track-result: non-finite numbers in ts_b", !flow.ts_b.is_finite());
    }

    pub fn pre_sync(&self, rough_delay: f64, ts_from: i64, ts_to: i64, search_step: f64, search_radius: f64) -> Option<(f64, f64)> {
        if self.problem.quats.is_empty() || search_step <= 0.0 || search_radius <= 0.0 {
            log::error!("Invalid params! quats.is_empty: {}, search_step: {search_step}, search_radius: {search_radius}", self.problem.quats.is_empty());
            return None;
        }
        let mut results = Vec::new();

        let timestamps: Vec<i64> = self.problem.frame_data.range(ts_from..ts_to).map(|(k, _)| *k).collect();

        let delays_len = (search_radius * 2.0) / search_step;
        let mut counter = 0.0;

        let mut delay = rough_delay - search_radius;
        while delay < rough_delay + search_radius {
            let cost: f64 = timestamps.par_iter().map(|ts| {
                let p = opt_compute_problem(*ts, delay, &self.problem);
                let m = opt_guess_translational_motion(&p, 20);
                let k = clamp_k(1.0 / (&p * &m).norm() * 1e2);
                let r = (&p * &m) * (k / m.norm());
                let rho = r.map(|v| libm::log1p(v * v).sqrt());
                rho.sum().sqrt()
            }).sum();
            results.push((cost, delay));
            delay += search_step;

            if let Some(ref cb) = self.progress_cb {
                if !cb((counter / delays_len) * 0.5) {
                    return None;
                }
            }
            counter += 1.0;
        }
        results.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Some(results[0])
    }

    pub fn sync(&self, initial_delay: f64, ts_from: i64, ts_to: i64, search_center: f64, search_radius: f64) -> Option<(f64, f64)> {
        if self.problem.quats.is_empty() {
            log::error!("Empty quats!");
            return None;
        }
        let mut gyro_delay = initial_delay;

        let costs: Vec<FrameState> = self.problem.frame_data.range(ts_from..ts_to).map(|(ts, _)| {
            let mut cost = FrameState::new(*ts, &self.problem);
            cost.motion_vec = cost.guess_motion(gyro_delay);
            cost.var_k = cost.guess_k(gyro_delay);
            cost
        }).collect();
        let costs = RefCell::new(costs);

        let mut delay_optimizer = Backtrack::default();
        delay_optimizer.set_hyper(2e-4, 0.1, 1e-3, 10);

        delay_optimizer.set_objective(|x: f64| -> (f64, f64) {
            costs.borrow().par_iter().map(|fs| {
                let mut cur_delay_g = 0.0;
                let mut tmp = DVector::from_element(0, 0.0);
                let cur_cost = fs.loss(x, &fs.motion_vec, &mut cur_delay_g, &mut tmp);
                (cur_cost, cur_delay_g)
            }).reduce_with(|a, b| (a.0 + b.0, a.1 + b.1)).unwrap_or_default()
        });

        let simple_objective = |x: f64| -> f64 {
            costs.borrow().par_iter().map(|fs| {
                fs.loss_single(x, &fs.motion_vec)
            }).sum()
        };

        delay_optimizer.set_objective_f_only(simple_objective);

        struct DelayOptInfo {
            step_size: f64
        }

        let delay_b = 0.3;
        let mut delay_v = 0.0;

        let mut converge_counter = 0;

        for _ in 0..400 {
            // Optimize motion
            {
                costs.borrow_mut().par_iter_mut().for_each(|fs| {
                    struct OptimizedFunction<'a> {
                        fs: &'a FrameState<'a>,
                        gyro_delay: f64,
                    }

                    impl<'a> Gradient for OptimizedFunction<'a> {
                        type Param = Vec<f64>;
                        type Gradient = Vec<f64>;

                        fn gradient(&self, w: &Self::Param) -> Result<Self::Param, Error> {
                            if w.iter().any(|v| !v.is_finite()) {
                                return Err(Error::msg("non-finite param"));
                            }
                            let mut del_jac = 0.0;
                            let mut grad = DVector::from_element(0, 0.0);
                            let _cost = self.fs.loss(self.gyro_delay, &DVector::from_column_slice(w), &mut del_jac, &mut grad);
                            Ok(grad.as_slice().to_vec())
                        }
                    }
                    impl<'a> CostFunction for OptimizedFunction<'a> {
                        type Param = Vec<f64>;
                        type Output = f64;

                        fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
                            if x.iter().any(|v| !v.is_finite()) {
                                return Err(Error::msg("non-finite param"));
                            }
                            Ok(self.fs.loss_single(self.gyro_delay, &DVector::from_column_slice(x)))
                        }
                    }

                    let cost = OptimizedFunction { fs: &fs, gyro_delay };
                    let linesearch = BacktrackingLineSearch::new(ArmijoCondition::new(1e-4).unwrap()).rho(0.5).unwrap();

                    let solver = LBFGS::new(linesearch, 10)
                        .with_tolerance_grad(1e-4).unwrap();

                    let executor = Executor::new(cost, solver)
                        .configure(|state| state.param(fs.motion_vec.as_slice().to_vec()).max_iters(200));

                    match executor.run() {
                        Ok(res) => {
                            // dbg!(&res.state().best_param);
                            // dbg!(&res.state().best_cost);
                            if let Some(ref best_param) = res.state().best_param {
                                fs.motion_vec = DVector::from_column_slice(&best_param);
                            } else {
                                log::error!("LBFGS error: bast_param is None");
                            }
                        },
                        Err(e) => {
                            log::error!("LBFGS error: {:?}", e);
                        }
                    }
                });
            }

            // Optimize delay
            let info = {
                let step = delay_optimizer.step(gyro_delay - delay_b * delay_v);

                delay_v = delay_b * delay_v + step;
                gyro_delay += delay_v;

                use argmin_math::ArgminL2Norm;

                DelayOptInfo { step_size: step.l2_norm() }
            };

            if info.step_size < 1e-4 {
                converge_counter += 1;
            } else {
                converge_counter = 0;
            }

            if converge_counter > 5 {
                break;
            }

            if (gyro_delay - search_center).abs() > search_radius {
                break;
            }
            // eprintln!("{} {}", gyro_delay, info.step_size);
        }

        Some((simple_objective(gyro_delay), gyro_delay))
    }

    pub fn full_sync(&self, initial_delay: f64, ts_from: i64, ts_to: i64, search_step: f64, search_radius: f64, iterations: usize) -> Option<(f64, f64)> {
        if let Some(mut delay) = self.pre_sync(initial_delay, ts_from, ts_to, search_step, search_radius) {
            for i in 0..iterations {
                if let Some(ref cb) = self.progress_cb {
                    if !cb(0.5 + (i as f64 / iterations as f64) * 0.5) {
                        return None;
                    }
                }
                if let Some(d) = self.sync(delay.1, ts_from, ts_to, initial_delay, search_radius) {
                    delay = d;
                }
            }
            if let Some(ref cb) = self.progress_cb {
                cb(1.0);
            }
            Some(delay)
        } else {
            None
        }
    }

    pub fn debug_pre_sync(&self, initial_delay: f64, ts_from: i64, ts_to: i64, search_radius: f64, delays: &mut [f64], costs: &mut [f64], point_count: usize) {
        let timestamps: Vec<i64> = self.problem.frame_data.range(ts_from..ts_to).map(|(k, _)| *k).collect();

        for i in 0..point_count {
            let delay = initial_delay - search_radius + 2.0 * search_radius * i as f64 / (point_count as f64 - 1.0);

            let cost: f64 = timestamps.par_iter().map(|ts| {
                let p = opt_compute_problem(*ts, delay, &self.problem);
                let m = opt_guess_translational_motion(&p, 20);
                let k = clamp_k(1.0 / (&p * &m).norm() * 1e2);
                let r = (&p * &m) * (k / m.norm());
                let rho = r.map(|v| libm::log1p(v * v).sqrt());
                rho.sum().sqrt()
            }).sum();
            delays[i] = delay;
            costs[i] = cost;
        }
    }
}


pub fn opt_compute_problem(timestamp_us: i64, gyro_delay: f64, data: &OptData) -> MatrixXx3<f64> {
    if let Some(flow) = data.frame_data.get(&timestamp_us) {
        let ap = &flow.rays_a;
        let bp = &flow.rays_b;
        let at = flow.ts_a.map(|v| (v - data.quats_start + gyro_delay) * data.sample_rate);
        let bt = flow.ts_b.map(|v| (v - data.quats_start + gyro_delay) * data.sample_rate);

        let mut problem = MatrixXx3::from_element(at.nrows(), 0.0);
        for i in 0..at.nrows() {
            let a = data.quats.eval(at[i]).normalize();
            let b = data.quats.eval(bt[i]).normalize();
            let ar = quat_rotate_point(&quat_conj(a), &ap[i]);
            let br = quat_rotate_point(&quat_conj(b), &bp[i]);

            problem.set_row(i, &ar.cross(&br).transpose());
        }

        problem
    } else {
        println!("frame not found {timestamp_us}: keys: {:?}", data.frame_data.keys());
        MatrixXx3::from_element(1, 0.0)
    }
}

pub fn opt_guess_translational_motion(problem: &MatrixXx3<f64>, max_iters: i32) -> DVector<f64> {
    if problem.is_empty() { return DVector::from_element(0, 0.0); }
    let mut nproblem = problem.clone();

    // TODO optimize
    for i in 0..nproblem.nrows() {
        let n = safe_normalize(Vector3::new(nproblem[(i, 0)], nproblem[(i, 1)], nproblem[(i, 2)]));
        nproblem[(i, 0)] = n[0];
        nproblem[(i, 1)] = n[1];
        nproblem[(i, 2)] = n[2];
    }

    let mut best_sol = DVector::from_element(3, 0.0);
    let mut least_med = f64::INFINITY;
    for _ in 0..max_iters {
        let rnd = mtrand(0, problem.nrows() - 1);
        let mut vs = [rnd, rnd];
        while problem.nrows() > 1 && vs[0] == vs[1] {
            vs[1] = mtrand(0, problem.nrows() - 1);
        }

        let v3 = problem.row(vs[0]).cross(&problem.row(vs[1]));
        let v: DVector<f64> = convert(safe_normalize(v3.transpose()));

        let residuals = &nproblem * &v;
        let mut residuals2 = residuals.map(|v| v * v).as_slice().to_vec(); // TODO

        residuals2.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let med = residuals2[residuals2.len() / 4];
        if med < least_med {
            least_med = med;
            best_sol = v;
        }
    }
    best_sol
}
