#include <iostream>
#include <ceres/ceres.h>
#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <fstream>
#include <chrono>
#include <assert.h>

class quat{
public:
    quat(double t, double x, double y, double z, double w):
    _t(t), _q(Eigen::Quaterniond(w, x, y, z)){}
    quat(double t, Eigen::Quaterniond q):
    _t(t), _q(q){}

    Eigen::Quaterniond _q;
    double _t;
};

class gyr{
public:
    gyr(double t, Eigen::Vector3d w):
    _t(t), _w(w){}
    double _t;
    Eigen::Vector3d _w;
};


Eigen::Quaterniond expmap(const Eigen::Vector3d &w) {
    Eigen::AngleAxisd aa(w.norm(), w.stableNormalized());
    Eigen::Quaterniond q;
    q = aa;
    return q;
}


std::vector<quat> loadData(const std::string& filename){
    std::fstream file(filename);
    std::vector<quat> qs;
    double t, x, y, z, w, _;
    while(file >> t >> _ >> _ >> _ >> x >> y >> z >> w){
        qs.emplace_back(quat(t, x, y, z ,w));
    }
    return qs;
}

std::vector<quat> loadGyrData(const std::string& filename){
    std::fstream file(filename);
    std::vector<quat> qs;
    std::vector<gyr> gs;
    double t, x, y, z;
    char c;
    while(file >> t >> c >> x >> c >> y >> c >> z){
        gs.emplace_back(gyr(t, Eigen::Vector3d(x, y, z)));
    }

    Eigen::Quaterniond cur_q(1, 0, 0, 0);
    for(size_t i = 0; i < gs.size() - 1; ++i){
        qs.emplace_back(quat(gs[i]._t, cur_q.x(), cur_q.y(), cur_q.z(), cur_q.w()));

        auto& g1 = gs[i];
        auto& g2 = gs[i + 1];
        double dt = g2._t - g1._t;
        Eigen::Vector3d r = (g1._w + g2._w) * dt / 2.0;
        cur_q = expmap(r) * cur_q;
    }
    return qs;
}



std::vector<Eigen::Quaterniond> interplate(const std::vector<quat>& q, const std::vector<double>& ts){
    std::vector<Eigen::Quaterniond> res;

    size_t pos = 0;
    for(size_t i = 0; i < ts.size(); ++i){
        while(pos < q.size() && ts[i] >= q[pos]._t) ++ pos;
        if(pos >= q.size()) break;
        auto q_start  = q[pos - 1];
        auto q_end    = q[pos];
        double ratio = (ts[i] - q_start._t) / (q_end._t - q_start._t);
        auto q_slerp = q_start._q.slerp(ratio, q_end._q);
        res.push_back(q_slerp);
    }

    return res;
}

void saveRotation(std::string path, const std::vector<double> logs){
    std::ofstream file(path);
    for(size_t i = 0; i < logs.size(); ++i){
        file << logs[i] << std::endl;
    }
    file.close();
}

void saveResult(const std::string& path, 
                const std::vector<double>& offset, 
                const std::vector<double>& tau,
                const std::vector<double>& ncc){
    std::ofstream file(path);
    for(size_t i = 0; i < offset.size(); ++i){
        file << offset[i] << " " << tau[i] << " " << ncc[i] << std::endl;
    }
    file.close();
}

double compute_NCC(const std::vector<quat>& qv, 
                    const std::vector<quat>& qi, 
                    double offset=0.0,
                    double init_offset=0.0,
                    double tau=0.5,
                    size_t patten_len=1000,
                    double start_t=1.0, 
                    double epsilon=0.1,
                    std::string save_name=""){

    double end_t = std::min(qi.back()._t - init_offset + offset, qv.back()._t) - tau - epsilon;

    std::vector<double> base_time;
    for(size_t i = 0; i < qi.size(); ++i){
        if(qi[i]._t - init_offset < start_t) continue;
        if(qi[i]._t - init_offset >= end_t) break;
        base_time.push_back(qi[i]._t);
        if(base_time.size() >= patten_len) break;
    }

    std::vector<double> vqvT1, vqvT2, vqiT1, vqiT2;
    for(auto& t: base_time) vqiT1.push_back(t);
    for(auto& t: base_time) vqiT2.push_back(t + tau);
    for(auto& t: base_time) vqvT1.push_back(t - init_offset + offset);
    for(auto& t: base_time) vqvT2.push_back(t - init_offset + offset + tau);

    auto qi1 = interplate(qi, vqiT1);
    auto qi2 = interplate(qi, vqiT2);
    auto qv1 = interplate(qv, vqvT1);
    auto qv2 = interplate(qv, vqvT2);
    
    std::vector<double> logs_qi, logs_qv;

    for(size_t i = 0; i < base_time.size(); ++i){

        Eigen::Quaterniond delta_qi = qi1[i].inverse() * qi2[i];
        Eigen::Quaterniond delta_qv = qv1[i].inverse() * qv2[i];
        Eigen::AngleAxisd log_qi(delta_qi);
        Eigen::AngleAxisd log_qv(delta_qv);

        logs_qi.push_back(log_qi.angle());
        logs_qv.push_back(log_qv.angle());
    }

    if (save_name != ""){
        saveRotation("./tmp/vicon_rotation_" + save_name + ".txt", logs_qv);
        saveRotation("./tmp/phone_rotation_" + save_name + ".txt", logs_qi);
    }

    double num1 ,num2 ,num3;
    num1 = num2 = num3 = 0.0;
    for(size_t i = 0; i < logs_qi.size(); i++){
        num1 += std::abs(logs_qi[i] * logs_qv[i]);
        num2 += std::abs(logs_qi[i] * logs_qi[i]);
        num3 += std::abs(logs_qv[i] * logs_qv[i]);
    }
    
    double ncc = num1 / std::sqrt(num2 * num3);
    return ncc;
}



int main(int argc, char** argv){

    if(argc < 3){
        std::cerr << "input: ./align_offset victo_gt_path estimated_traj_path" << std::endl;
        exit(0);
    }
    std::string qv_path = argv[1];
    std::string qi_path = argv[2];

    std::cout << std::endl;
    std::cout << "--vicon traj: " << qv_path << std::endl;
    std::cout << "--phone traj: " << qi_path << std::endl;
    std::cout << std::endl;

    std::vector<quat> qv = loadData(qv_path);
    std::vector<quat> qi = loadData(qi_path);
    // std::vector<quat> qi = loadGyrData(qi_path);


    std::cout << "--vicon traj data size: " << qv.size() << std::endl;
    std::cout << "--phone traj data size: " << qi.size() << std::endl;
    
    double init_offset = qi[0]._t - qv[0]._t;
    double offset_lim[2] = {-1.0, 5.0};
    double taus[5] = {1, 2, 3, 4, 5};

    std::vector<double> vOffset, vNCC, vTau;

    double max_ncc = 0;
    double best_offset = -10;
    for(double tau = 1; tau < 6; ++tau){
        for(double offset = offset_lim[0]; offset < offset_lim[1]; offset += 0.001){
            double ncc = compute_NCC(qv, qi, offset, init_offset, tau);
            printf("-- offset: %f, tau: %f, NCC: %f\n", offset, tau, ncc);
            vOffset.push_back(offset);
            vNCC.push_back(ncc);
            vTau.push_back(tau);
            if (tau == 1 && max_ncc < ncc){
                max_ncc = ncc;
                best_offset = offset;
            }
        }
    }

    std::cout << std::endl;
    std::cout << "--max ncc         :" << max_ncc << std::endl;
    std::cout << "--best offset     :" << best_offset << std::endl;


    compute_NCC(qv, qi, 0, init_offset, 1.0, 1000, 1.0, 0.1, "before");

    double set_offet = best_offset;
    compute_NCC(qv, qi, set_offet, init_offset, 1.0, 1000, 1.0, 0.1, "after");

    std::string save_dir = "./tmp/result_offset_ncc.txt";
    saveResult(save_dir, vOffset, vTau, vNCC);

    return 0;
}