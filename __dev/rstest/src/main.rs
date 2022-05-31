use rs_sync::SyncProblem;
mod cpp_wrapper;

fn main() {
    let point_count = 200;
    let mut delays = Vec::new();
    let mut costs = Vec::new();
    delays.resize(point_count, 0.0);
    costs.resize(point_count, 0.0);

    walkdir::WalkDir::new("test_data/").into_iter().flatten().for_each(|entry| {
        let path = entry.path().to_string_lossy();
        if !path.ends_with(".bin") { return; }
        
        println!("{path}:");
        let d = cpp_wrapper::load_data_from_file(&path);

        { // C++ wrapper
            let mut sync = cpp_wrapper::SyncProblem::new();

            sync.set_gyro_quaternions(&d.timestamps, &d.quats);
            for x in &d.perframe {        
                sync.set_track_result(x.timestamp_us, &x.tsa, &x.tsb, &x.pointsa, &x.pointsb);
            }
            sync.debug_pre_sync(d.initial_delay / 1000.0, d.from_ts, d.to_ts, d.presync_radius / 1000.0, delays.as_mut_slice(), costs.as_mut_slice(), point_count);

            let mut delay = sync.pre_sync(d.initial_delay / 1000.0, d.from_ts, d.to_ts, d.presync_step / 1000.0, d.presync_radius / 1000.0);
            println!("cpp pre_sync: {:.4}", delay.1 * 1000.0);
            for _ in 0..4 {
                delay = sync.sync(delay.1, d.from_ts, d.to_ts);
            }
            println!("cpp offset: {:.4}", delay.1 * 1000.0);
        }

        { // Rust native
            let mut sync = SyncProblem::new();

            sync.set_gyro_quaternions(&d.timestamps, &d.quats);
            for x in &d.perframe {        
                sync.set_track_result(x.timestamp_us, &x.tsa, &x.tsb, &x.pointsa, &x.pointsb);
            }
            sync.debug_pre_sync(d.initial_delay / 1000.0, d.from_ts, d.to_ts, d.presync_radius / 1000.0, delays.as_mut_slice(), costs.as_mut_slice(), point_count);

            let mut delay = sync.pre_sync(d.initial_delay / 1000.0, d.from_ts, d.to_ts, d.presync_step / 1000.0, d.presync_radius / 1000.0);
            println!("rust pre_sync: {:.4}", delay.1 * 1000.0);
            for _ in 0..4 {
                delay = sync.sync(delay.1, d.from_ts, d.to_ts);
            }
            println!("rust offset: {:.4}", delay.1 * 1000.0);
        }

        let mut out = String::new();
        for i in 0..point_count {
            out.push_str(&format!("{},{}\n", delays[i], costs[i]))
        }
        std::fs::write("dbg.csv", out).unwrap();
    });

}

