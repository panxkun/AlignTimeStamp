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

std::vector<quat> loadData(const std::string& filename){
    std::fstream file(filename);
    std::vector<quat> qs;
    float t, x, y, z, w, _;
    while(file >> t >> _ >> _ >> _ >> x >> y >> z >> w){
        qs.emplace_back(quat(t, x, y, z ,w));
    }
    
    return qs;
}

std::vector<double> process_timestamps(std::vector<quat>& qs, double start=5.0){
    std::vector<double> res;
    double current_time = -1.0;
    for(auto& q: qs){
        if(current_time < 0){
            current_time = q._t;
            continue;
        }
        else if(q._t - current_time < start){
            continue;
        }else{
            res.push_back(q);
        }
    }
    return res;
}

Eigen::Quaterniond Quaternion_S_lerp(Eigen::Quaterniond start_q, Eigen::Quaterniond end_q, double t)
{
    Eigen::Quaterniond lerp_q;
    double cos_angle = start_q.x() * end_q.x() + start_q.y() * end_q.y() + start_q.z() * end_q.z() + start_q.w() * end_q.w();

    // If the dot product is negative, the quaternions have opposite handed-ness and slerp won't take
    // the shorter path. Fix by reversing one quaternion.
    if (cos_angle < 0) {
        end_q.x() = -end_q.x();
        end_q.y() = -end_q.y();
        end_q.z() = -end_q.z();
        end_q.w() = -end_q.w();
        cos_angle = -cos_angle;
    }

    double ratio_A, ratio_B;
    if (cos_angle > 0.99995f) {
        ratio_A = 1.0f - t;
        ratio_B = t;
    }
    else {
        double sin_angle = sqrt( 1.0f - cos_angle * cos_angle);
        double angle = atan2(sin_angle, cos_angle);
        ratio_A = sin((1.0f - t) * angle)  / sin_angle;
        ratio_B = sin(t * angle) / sin_angle;
    }

    lerp_q.x() = ratio_A * start_q.x() + ratio_B * end_q.x();
    lerp_q.y() = ratio_A * start_q.y() + ratio_B * end_q.y();
    lerp_q.z() = ratio_A * start_q.z() + ratio_B * end_q.z();
    lerp_q.w() = ratio_A * start_q.w() + ratio_B * end_q.w();

    return lerp_q.normalized();
}


std::vector<quat> interpolate(std::vector<quat> q, std::vector<double> ts){
    
    assert(ts > 0);

    std::vector<quat> res;

    size_t pos = 0;
    size_t len = ts.size();
    for(size_t i = 0; i < len; ++i){
        double& cur_t = ts[i];
        if(cur_t < q[0]._t){
            res.push_back(quat(cur_t, q[0]._q));
            continue;
        }
        else if(cur_t > q[len - 1]._t){

        }

        while(pos < q.size() && ts[i] > q[pos]._t) ++ pos;
        if(pos >= q.size()) break;
        auto q_start  = q[pos];
        while(pos < q.size() && ts[i] < q[pos]._t) ++ pos;
        if(pos >= q.size()) break;
        auto q_end  = q[pos];

        std::cout << q_start._t << " " << ts[i] << " " << q_end._t << std::endl;
        double ratio = (ts[i] - q_start._t) / (q_end._t - q_start._t);
        auto q_slerp = Quaternion_S_lerp(q_start._q, q_end._q, ratio);
        res.push_back(quat(ts[i], q_slerp));
    }

    std::cout << "aftre slerp: " << res.size() << std::endl;

    return res;
}

struct RelativeRoatationError{
    RelativeRoatationError(const std::vector<quat>& qv, const std::vector<quat>& qi, double tau):
    _qv(qv), _qi(qi), _tau(tau){
    }

    // template<typename T>
    bool operator() (const double* const t, double* residual) const  
    {
        std::vector<double> vqvT1, vqvT2, vqiT1, vqiT2;
        std::cout << "est t: " << *t << std::endl;

        for(auto& q: _qi) vqiT1.push_back(q._t);
        for(auto& q: _qi) vqiT2.push_back(q._t + _tau);
        for(auto& q: _qv) vqvT1.push_back(t[0] + q._t);
        for(auto& q: _qv) vqvT2.push_back(t[0] + q._t + _tau);

        std::vector<quat> qi_t1, qv_t1;
        std::vector<quat> qi_t1 = interpolate(_qi, vqiT1);
        std::vector<quat> qi_t2 = interpolate(_qi, vqiT2);
        std::vector<quat> qv_t1 = interpolate(_qv, vqvT1);
        std::vector<quat> qv_t2 = interpolate(_qv, vqvT2);

        double error = 0.0;

        for(size_t i = 0; i < vqvT1.size(); ++i){
            Eigen::Quaterniond delta_qv = qv_t1[i]._q.inverse() * qv_t2[i]._q;
            Eigen::Quaterniond delta_qi = qi_t1[i]._q.inverse() * qi_t2[i]._q;
            Eigen::Quaterniond log_res = delta_qv.inverse() * delta_qi;
            Eigen::AngleAxisd aa(log_res);
            error += aa.angle() * aa.angle();
        }

        residual[0] = error;

        return true;
    }

    static ceres::CostFunction* Create(const std::vector<quat>& qv, const std::vector<quat>& qi, double tau){
        return (new ceres::NumericDiffCostFunction<RelativeRoatationError, ceres::CENTRAL, 1, 1>(
                new RelativeRoatationError(qv, qi, tau)));
   }

    std::vector<quat> _qv;
    std::vector<quat> _qi;
    double _tau;
};

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

    std::vector<quat> qv = loadData(qv_path);
    std::vector<quat> qi = loadData(qi_path);

    std::cout << "--vicon traj data size: " << qv.size() << std::endl;
    std::cout << "--phone traj data size: " << qi.size() << std::endl;
    
    qi = process_timestamps(qi);
    double t[1] = {qi[0]._t - qt[0]._t};
    double tau = 0.1;

    google::InitGoogleLogging(argv[0]);

    ceres::Problem problem;
    ceres::CostFunction* cost_function = RelativeRoatationError::Create(qv, qi, tau);
    problem.AddResidualBlock(cost_function, nullptr, t);
    
    ceres::Solver::Options options;    
    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;  
    options.minimizer_progress_to_stdout = true;  

    ceres::Solver::Summary summary;                
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    ceres::Solve ( options, &problem, &summary);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>( t2-t1 );
    std::cout << "solve time cost = " << time_used.count() << " seconds. " << std::endl;

    return 0;
}