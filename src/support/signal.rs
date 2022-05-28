// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright © 2021-2022 Vladimir Pinchuk (https://github.com/VladimirP1)
// Ported to Rust by Adrian <adrian.eddy at gmail>

use nalgebra::*;
/*
pub fn gyro_lowpass(samples: &mut Matrix, divider: i32) {
    if divider < 2 { return; }
    let ita = 1.0 / (std::f64::consts::PI / divider as f64).tan();
    let q = 2.0.sqrt();
    let b0 = 1.0 / (1.0 + q * ita + ita * ita);
    let b1 = 2 * b0;
    let b2 = b0;
    let a1 = 2.0 * (ita * ita - 1.0) * b0;
    let a2 = -(1.0 - q * ita + ita * ita) * b0;

    let mut out = [samples.col(0), samples.col(1), samples.col(2)];
    for i in 2..samples.ncols() {
        out[2] = b0 * samples.col(i) + b1 * samples.col(i - 1) + b2 * samples.col(i - 2) +
                 a1 * out[2 - 1] + a2 * out[2 - 2];
        samples.col(i - 2) = out[0];
        // left shift
        out[0] = out[1];
        out[1] = out[2];
    }
    // reverse pass
    out[0] = samples.col(samples.ncols() - 1);
    out[1] = samples.col(samples.ncols() - 2);
    for j in 2..samples.ncols() {
        let i = samples.ncols() - j - 1;
        out[2] = b0 * samples.col(i) + b1 * samples.col(i + 1) + b2 * samples.col(i + 2) +
                 a1 * out[2 - 1] + a2 * out[2 - 2];
        samples.col(i + 2) = out[0];
        // left shift
        out[0] = out[1];
        out[1] = out[2];
    }
}

pub fn gyro_upsample(samples: &mut Matrix, multiplier: usize) {
    if multiplier < 2 { return; }
    let length_new = samples.ncols() * multiplier;
    let length = length_new / multiplier;
    let half_mult = multiplier / 2;
    let old_samples_base = length_new - length;
    // std::copy_n(samples, length, samples + old_samples_base);
    samples.cols(old_samples_base, old_samples_base + length - 1) = samples.cols(0, length - 1);

    for i in 0..length_new {
        if (i + half_mult) % multiplier {
            samples.col(i) = {};
        } else {
            samples.col(i) = samples.col(i / multiplier + old_samples_base);
        }
    }

    gyro_lowpass(samples, multiplier * 4);
}

pub fn gyro_decimate(samples: &mut Matrix, divider: usize) {
    if divider < 2 { return; }
    let samples2 = samples.clone();
    samples.resize(samples.nrows(), samples.ncols() / divider);
    for i in 0..samples.ncols() / divider {
        samples.col(i) = samples2.col(i * divider);
    }
}

pub fn gyro_interpolate(timestamps: &mut Matrix, gyro: &mut Matrix) -> i32 {
    let actual_sr = timestamps.len() / (timestamps.back() - timestamps.front());
    let rounded_sr = ((actual_sr / 50.0).round() * 50) as i32;

    std::vector<double> new_timestamps_vec;
    for (double sample = std::ceil(timestamps.front() * rounded_sr);
         sample / rounded_sr < timestamps.back(); sample += 1)
        new_timestamps_vec.push_back(sample / rounded_sr);

    arma::mat new_timestamps(new_timestamps_vec.data(), 1, new_timestamps_vec.size());
    arma::mat new_gyro(3, new_timestamps_vec.size());
    arma::mat tmp;
    arma::interp1(timestamps, gyro.row(0), new_timestamps, tmp);
    new_gyro.row(0) = tmp;
    arma::interp1(timestamps, gyro.row(1), new_timestamps, tmp);
    new_gyro.row(1) = tmp;
    arma::interp1(timestamps, gyro.row(2), new_timestamps, tmp);
    new_gyro.row(2) = tmp;

    gyro = new_gyro;
    timestamps = new_timestamps;

    rounded_sr
}
*/