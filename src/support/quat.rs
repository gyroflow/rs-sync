// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright © 2021-2022 Vladimir Pinchuk (https://github.com/VladimirP1)
// Ported to Rust by Adrian <adrian.eddy at gmail>

#![allow(unused)]

use nalgebra::*;

pub fn quat_from_aa(aa: &Vector3<f64>) -> Vector4<f64> {
    let theta_squared = aa.dot(aa);

    if theta_squared > 0.0 {
        let theta = theta_squared.sqrt();
        let half_theta = theta * 0.5;
        let k = half_theta.sin() / theta;
        Vector4::new(half_theta.cos(), aa[0] * k, aa[1] * k, aa[2] * k)
    } else {
        let k = 0.5;
        Vector4::new(1.0, aa[0] * k, aa[1] * k, aa[2] * k)
    }
}

pub fn quat_to_aa(q: &Vector4<f64>) -> Vector3<f64> {
    let xyz = q.fixed_rows::<3>(1);
    let sin_squared_theta = xyz.dot(&xyz);

    if sin_squared_theta <= 0.0 {
        return xyz * 2.0;
    }

    let sin_theta = sin_squared_theta.sqrt();
    let cos_theta = q[0];
    let two_theta = 2.0 * if cos_theta < 0.0 { (-sin_theta).atan2(-cos_theta) } else { sin_theta.atan2(cos_theta) };
    let k = two_theta / sin_theta;
    xyz * k
}

pub fn quat_prod(p: &Vector4<f64>, q: &Vector4<f64>) -> Vector4<f64> {
    Vector4::new(p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
                 p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
                 p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
                 p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0])
}

pub fn quat_conj(mut q: Vector4<f64>) -> Vector4<f64> {
    *q.fixed_rows_mut::<3>(1) = *-q.fixed_rows::<3>(1);
    q
}

pub fn quat_rotate_point(q: &Vector4<f64>, p: &Vector3<f64>) -> Vector3<f64> {
    Vector3::from(quat_prod(q, &quat_prod(&Vector4::new(0.0, p[0], p[1], p[2]), &quat_conj(*q))).fixed_rows::<3>(1))
}

pub fn quat_double(p: &Vector4<f64>, q: &Vector4<f64>) -> Vector4<f64> {
    2.0 * p.dot(q) * q - p
}

pub fn quat_bisect(p: &Vector4<f64>, q: &Vector4<f64>) -> Vector4<f64> { (p + q) * 0.5 }

pub fn quat_slerp(p: &Vector4<f64>, mut q: Vector4<f64>, t: f64) -> Vector4<f64> {
    if p.dot(&q) < 0.0 {
        q = -q;
    }

    let theta = p.dot(&q).acos();

    // TODO: check if differentiable
    let (mult1, mult2) = if theta > 1e-9 {
        let sin_theta = theta.sin();
        (
            ((1.0 - t) * theta).sin() / sin_theta,
            (t * theta).sin() / sin_theta
        )
    } else {
        (1.0 - t, t)
    };

    mult1 * p + mult2 * q
}

pub fn quat_squad(p0: &Vector4<f64>, p1: &Vector4<f64>, p2: &Vector4<f64>, p3: &Vector4<f64>, t: f64) -> Vector4<f64> {
    let a0 = quat_bisect(&quat_double(p0, p1), p2);
    let a1 = quat_bisect(&quat_double(p1, p2), p3);
    let b1 = quat_double(&a1, p2);
    let i0 = p1;
    let mut i1 = a0;
    let mut i2 = b1;
    let i3 = p2;
    i1 = (i1 + 2.0 * i0) / 3.0;
    i2 = (i2 + 2.0 * i3) / 3.0;
    let j0 = quat_slerp(i0, i1, t);
    let j1 = quat_slerp(&i1, i2, t);
    let j2 = quat_slerp(&i2, *i3, t);
    quat_slerp(&quat_slerp(&j0, j1, t), quat_slerp(&j1, j2, t), t)
}

pub fn quat_lerp(p: &Vector4<f64>, q: &Vector4<f64>, t: f64) -> Vector4<f64> { p * (1.0 - t) + q * t }

pub fn quat_quad(p0: &Vector4<f64>, p1: &Vector4<f64>, p2: &Vector4<f64>, p3: &Vector4<f64>, t: f64) -> Vector4<f64> {
    let mut a0 = quat_bisect(&quat_double(p0, p1), p2);
    let     a1 = quat_bisect(&quat_double(p1, p2), p3);
    let mut b1 = quat_double(&a1, p2);
    a0 = (a0 + 2.0 * p1) / 3.0;
    b1 = (b1 + 2.0 * p2) / 3.0;
    let j0 = quat_lerp(p1, &a0, t);
    let j1 = quat_lerp(&a0, &b1, t);
    let j2 = quat_lerp(&b1, p2, t);
    quat_lerp(&quat_lerp(&j0, &j1, t), &quat_lerp(&j1, &j2, t), t)
}
