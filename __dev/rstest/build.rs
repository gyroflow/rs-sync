fn main() {
    let mut config = cpp_build::Config::new();

    let path = format!("{}/../cpp/rs-sync", std::env::var("CARGO_MANIFEST_DIR").unwrap());

    if std::env::var("CARGO_CFG_TARGET_ENV").unwrap() == "msvc" {
        config.flag("/Zc:__cplusplus");
        config.flag("/std:c++17");
    } else {
        config.flag("-std=c++17");
    }
    config.define("_USE_MATH_DEFINES", None);
    config.define("ARMA_DONT_USE_BLAS", None);
    config.define("ARMA_DONT_USE_LAPACK", None);

    config.include(&format!("{}/src/core/public", path));
    config.include(&format!("{}/src/core", path));
    config.include(&format!("{}/src/core_support", path));
    config.include(&format!("{}/../armadillo/include", path));
    config.include(&format!("{}/../ensmallen/include", path));

    config.file(&format!("{}/src/core_support/backtrack.cpp", path));
    config.file(&format!("{}/src/core_support/minispline.cpp", path)); 
    config.file(&format!("{}/src/core_support/ndspline.cpp", path));     
    config.file(&format!("{}/src/core_support/quat.cpp", path));
    config.file(&format!("{}/src/core_support/panic.cpp", path));
    config.file(&format!("{}/src/core/core_private.cpp", path));
       
    config.build("src/main.rs");

    println!("cargo:rerun-if-changed={}/src/", path);
    // println!("cargo:rustc-link-lib=tbb");
}
