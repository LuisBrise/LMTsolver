//**********************************************************************
// Generate paths
//********************************************************************** 

std::string get_current_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y%m%d_%H%M%S");
    return ss.str();
}

void write_metadata(const fsm::path& dir, const SimulationConfig& config) {
    std::ofstream meta(dir / "simulation_parameters.txt");
    meta << "Material: " << config.material << "\n"
         << "Lmax: " << config.Lmax << "\n"
         << "Velocity: " << config.velocity << "c\n"
         << "Integration surface: " << config.r << " nm\n"
         << "Impact param: " << config.b << " nm\n"
         << "NP radius: " << config.a << " nm\n"
         << "Date: " << get_current_timestamp();
}

fsm::path generate_multconvan_path(const SimulationConfig& config) {
    
    // Create main directory path
    fsm::path dir = "results";
    dir /= ("material=" + config.material);

    // Handle decimal precision for 'a'
    std::ostringstream oss_a;
    oss_a << "a=" << std::fixed << std::setprecision(2) << config.a << "nm";
    dir /= oss_a.str();  // e.g., "a=50.50nm"

    dir /= config.timestamp;

        // Parameter categorization
    if (config.isVScan) {
        std::ostringstream oss_b;
        oss_b << "v_scan_for_b=" << std::setprecision(2) << config.b << "nm";
        dir /= oss_b.str();
    } 
    else if (config.isBScan) {
        std::ostringstream oss_v;
        oss_v << "b_scan_for_v=" << std::setprecision(2) << config.velocity << "c";
        dir /= oss_v.str();
    }
    else if (config.isBvsVContour) {
        std::ostringstream oss_vb;
        oss_vb  << "b_vs_v_contour_for_v=" << std::setprecision(2) << config.velocity << "c"
           << "_and_b=" << std::setprecision(2) << config.b << "nm";
        dir /= oss_vb.str();
    }

    dir /= "MultipolarConvergence";  

    // Create directories
    fsm::create_directories(dir);

    // Generate filename
    std::ostringstream filename;
    filename << "MultipolarConvergence_"
             << config.material << "_"
             << "v" << std::fixed << std::setprecision(2) << config.velocity << "c_"
             << "b" << std::fixed << std::setprecision(2) << config.b << "nm_"
             << "a" << std::fixed << std::setprecision(2) << config.a << "nm_"
             << ".dat";

    write_metadata(dir, config);
    // Combine directory and filename
    return dir / filename.str();
}

fsm::path generate_AMT_path(const SimulationConfig& config,
                           const std::string& component = "",
                           const std::string& error = "")  {
    
    // Create main directory path
    fsm::path dir = "results";
    dir /= ("material=" + config.material);

    // Handle decimal precision for 'a'
    std::ostringstream oss_a;
    oss_a << "a=" << std::fixed << std::setprecision(2) << config.a << "nm";
    dir /= oss_a.str();  // e.g., "a=50.50nm"

    //Get timestamp
    dir /= config.timestamp;

    // Generate filename
    std::ostringstream filename;
    filename << error; // write const std::string& error if given, if not leaves empty

    // Parameter categorization
    if (config.isVScan) {
        std::ostringstream oss_b;
        oss_b << "v_scan_for_b=" << std::setprecision(2) << config.b << "nm";
        dir /= oss_b.str();
        filename << "DL" << component << "_"
                 << config.material << "_"
                 << "b" << std::setprecision(2) << config.b << "nm_"
                 << "a" << std::setprecision(2) << config.a << "nm_"
                 << "L" << config.Lmax << "_"
                 << ".dat";
    } 
    else if (config.isBScan) {
        std::ostringstream oss_v;
        oss_v << "b_scan_for_v=" << std::setprecision(2) << config.velocity << "c";
        dir /= oss_v.str();
        filename << "DL" << component << "_"
                 << config.material << "_"
                 << "v" << std::setprecision(2) << config.velocity << "c_"
                 << "a" << std::setprecision(2) << config.a << "nm_"
                 << "L" << config.Lmax << "_"
                 << ".dat";
    }
    else if (config.isBvsVContour) {
        std::ostringstream oss_vb;
        oss_vb  << "b_vs_v_contour_for_v=" << std::setprecision(2) << config.velocity << "c"
           << "_and_b=" << std::setprecision(2) << config.b << "nm";
        dir /= oss_vb.str();
        filename << "DL" << component << "_"
                 << config.material << "_"
                 << "a" << std::setprecision(2) << config.a << "nm_"
                 << "L" << config.Lmax << "_"
                 << ".dat";
    }
    else{
        filename << "DL" << component << "_"
                 << config.material << "_"
                 << "v" << std::setprecision(2) << config.velocity << "c_"
                 << "b" << std::setprecision(2) << config.b << "nm_"
                 << "a" << std::setprecision(2) << config.a << "nm_"
                 << "L" << config.Lmax << "_"
                 << ".dat";
    }
    // Lmax subdirectory
    dir /= ("Lmax=" + std::to_string(config.Lmax));

    // Create directories (including parents)
    fsm::create_directories(dir);
    write_metadata(dir, config);
    // Combine directory and filename
    return dir / filename.str();
}

fsm::path generate_dldw_path(const SimulationConfig& config,
                             const std::string& component = "")  {
    
    // Create main directory path
    fsm::path dir = "results";
    dir /= ("material=" + config.material);

    // Handle decimal precision for 'a'
    std::ostringstream oss_a;
    oss_a << "a=" << std::fixed << std::setprecision(2) << config.a << "nm";
    dir /= oss_a.str();  // e.g., "a=50.50nm"

    //timestamp
    dir /= config.timestamp;
    
    // Parameter categorization
    if (config.isVScan) {
        std::ostringstream oss_b;
        oss_b << "v_scan_for_b=" << std::setprecision(2) << config.b << "nm";
        dir /= oss_b.str();
    } 
    else if (config.isBScan) {
        std::ostringstream oss_v;
        oss_v << "b_scan_for_v=" << std::setprecision(2) << config.velocity << "c";
        dir /= oss_v.str();
    }
    else if (config.isBvsVContour) {
        std::ostringstream oss_vb;
        oss_vb  << "b_vs_v_contour_for_v=" << std::setprecision(2) << config.velocity << "c"
           << "_and_b=" << std::setprecision(2) << config.b << "nm";
        dir /= oss_vb.str();
    }

    // Lmax subdirectory
    dir /= ("Lmax=" + std::to_string(config.Lmax));
    
    // Generate filename
    std::ostringstream filename;
    filename << ("dldw" + component + "_")
             << config.material << "_"
             << "v" << std::setprecision(2) << config.velocity << "c_"
             << "b" << std::setprecision(2) << config.b << "nm_"
             << "a" << std::setprecision(2) << config.a << "nm_"
             << config.timestamp << ".dat";

    // Create directories (including parents)
    fsm::create_directories(dir);

    // Filename
    return dir / filename.str();
}