#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <spdlog/spdlog.h>
#include <iostream>
#include "Labs/4-Animation/tasks.h"
#include "IKSystem.h"
#include "CustomFunc.inl"


namespace VCX::Labs::Animation {
    void ForwardKinematics(IKSystem & ik, int StartIndex) {
        if (StartIndex == 0) {
            ik.JointGlobalRotation[0] = ik.JointLocalRotation[0];
            ik.JointGlobalPosition[0] = ik.JointLocalOffset[0];
            StartIndex                = 1;
        }
        
        for (int i = StartIndex; i < ik.JointLocalOffset.size(); i++) {
            // your code here: forward kinematics
            ik.JointGlobalPosition[i] = ik.JointGlobalPosition[i-1] + ik.JointGlobalRotation[i-1] * ik.JointLocalOffset[i];
            ik.JointGlobalRotation[i] = ik.JointGlobalRotation[i-1] * ik.JointLocalRotation[i];
        }
    }

    void InverseKinematicsCCD(IKSystem & ik, const glm::vec3 & EndPosition, int maxCCDIKIteration, float eps) {
        ForwardKinematics(ik, 0);
        // These functions will be useful: glm::normalize, glm::rotation, glm::quat * glm::quat
        for (int CCDIKIteration = 0; CCDIKIteration < maxCCDIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; CCDIKIteration++) {
            // your code here: ccd ik
            int njoints = ik.NumJoints();
            for (int i = njoints-2; i >= 0; --i){
                glm::vec3 v1 = ik.JointGlobalPosition[njoints-1] - ik.JointGlobalPosition[i];
                glm::vec3 v2 = EndPosition - ik.JointGlobalPosition[i];
                Eigen::Vector3d vv1 {v1.x, v1.y, v1.z};
                Eigen::Vector3d vv2 {v2.x, v2.y, v2.z};
                Eigen::Quaterniond a = Eigen::Quaterniond::FromTwoVectors(vv1, vv2);
                glm::quat delta (a.w(), a.x(), a.y(), a.z());
                ik.JointLocalRotation[i] = delta * ik.JointLocalRotation[i];
                ForwardKinematics(ik, i);
            }
        }
    }

    void InverseKinematicsFABR(IKSystem & ik, const glm::vec3 & EndPosition, int maxFABRIKIteration, float eps) {
        ForwardKinematics(ik, 0);
        int nJoints = ik.NumJoints();
        std::vector<glm::vec3> backward_positions(nJoints, glm::vec3(0, 0, 0)), forward_positions(nJoints, glm::vec3(0, 0, 0));
        for (int IKIteration = 0; IKIteration < maxFABRIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; IKIteration++) {
            // task: fabr ik
            // backward update
            glm::vec3 next_position         = EndPosition;
            backward_positions[nJoints - 1] = EndPosition;

            for (int i = nJoints - 2; i >= 0; i--) {
                glm::vec3 vec = ik.JointGlobalPosition[i] - next_position;
                float ll     = sqrt(glm::dot(vec, vec));
                next_position = next_position + (ik.JointOffsetLength[i + 1] * vec / ll);
                backward_positions[i] = next_position;
            }

            // forward update
            glm::vec3 now_position = ik.JointGlobalPosition[0];
            forward_positions[0] = ik.JointGlobalPosition[0];
            for (int i = 0; i < nJoints - 1; i++) {
                // your code here
                glm::vec3 vec = backward_positions[i + 1] - now_position;
                float    ll  = sqrt(glm::dot(vec, vec));
                now_position  = now_position + (ik.JointOffsetLength[i + 1] * vec / ll);
                forward_positions[i + 1] = now_position;
            }
            ik.JointGlobalPosition = forward_positions; // copy forward positions to joint_positions
        }

        // Compute joint rotation by position here.
        for (int i = 0; i < nJoints - 1; i++) {
            ik.JointGlobalRotation[i] = glm::rotation(glm::normalize(ik.JointLocalOffset[i + 1]), glm::normalize(ik.JointGlobalPosition[i + 1] - ik.JointGlobalPosition[i]));
        }
        ik.JointLocalRotation[0] = ik.JointGlobalRotation[0];
        for (int i = 1; i < nJoints - 1; i++) {
            ik.JointLocalRotation[i] = glm::inverse(ik.JointGlobalRotation[i - 1]) * ik.JointGlobalRotation[i];
        }
        ForwardKinematics(ik, 0);
    }

    IKSystem::Vec3ArrPtr IKSystem::BuildCustomTargetPosition() {
        // get function from https://www.wolframalpha.com/input/?i=Albert+Einstein+curve
        // int nums = 5000;
        // using Vec3Arr = std::vector<glm::vec3>;
        // std::shared_ptr<Vec3Arr> custom(new Vec3Arr(nums));
        // int index = 0;
        // for (int i = 0; i < nums; i++) {
        //     float theta = (2*3.14*i) / nums;
        //     float x_val = 16 * pow(sin(theta), 3) / 35 + 0.5f;
        //     float y_val = (13 * cos(theta) - 5*cos(2*theta) - 2*cos(3*theta) - cos(4*theta)) / 35 + 0.5f;
        //     if (std::abs(x_val) < 1e-3 || std::abs(y_val) < 1e-3) continue;
        //     (*custom)[index++] = glm::vec3(1.6f - x_val, 0.0f, y_val - 0.2f);
        // }
        // custom->resize(index);
        using Vec3Arr = std::vector<glm::vec3>;
        std::shared_ptr<Vec3Arr> custom(new Vec3Arr(1000));
        int                      index = 0;
        const int max_count = 100;
        for (int i = 0; i < max_count; ++i){
            float theta = (2 * 3.14 * i) / max_count;
            (*custom)[index++] = glm::vec3(16 * pow(sin(theta), 3) / 35 + 0.5f, 0, (13 * cos(theta) - 5*cos(2*theta) - 2*cos(3*theta) - cos(4*theta)) / 35 + 0.5f);
        }
        custom->resize(index);
        return custom;
    }

    void AdvanceMassSpringSystem(MassSpringSystem & system, float const dt) {
        // your code here: rewrite following code
        // int const steps = 1000;
        // float const ddt = dt / steps; 
        // for (std::size_t s = 0; s < steps; s++) {
        //     std::vector<glm::vec3> forces(system.Positions.size(), glm::vec3(0));
        //     for (auto const spring : system.Springs) {
        //         auto const p0 = spring.AdjIdx.first;
        //         auto const p1 = spring.AdjIdx.second;
        //         glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
        //         glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
        //         glm::vec3 const e01 = glm::normalize(x01);
        //         glm::vec3 f = (system.Stiffness * (glm::length(x01) - spring.RestLength) + system.Damping * glm::dot(v01, e01)) * e01;
        //         forces[p0] += f;
        //         forces[p1] -= f;
        //     }
        //     for (std::size_t i = 0; i < system.Positions.size(); i++) {
        //         if (system.Fixed[i]) continue;
        //         system.Velocities[i] += (glm::vec3(0, -system.Gravity, 0) + forces[i] / system.Mass) * ddt;
        //         system.Positions[i] += system.Velocities[i] * ddt;
        //     }
        // }
        int const steps = 1;
        float const ddt = dt / steps; 
        for (std::size_t s = 0; s < steps; s++) {
            std::vector<glm::vec3> forces(system.Positions.size(), glm::vec3(0));
            int n = system.Positions.size();
            Eigen::SparseMatrix<float> M(3 * n, 3 * n);
            std::vector<Eigen::Triplet<float>> coefficients;
            coefficients.clear();
            for (auto const spring : system.Springs) {
                auto const p0 = spring.AdjIdx.first;
                auto const p1 = spring.AdjIdx.second;
                glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
                glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
                glm::vec3 e01 = glm::normalize(x01);
                glm::vec3 f = (system.Stiffness * (glm::length(x01) - spring.RestLength) + system.Damping * glm::dot(v01, e01)) * e01;
                glm::vec3 vec = -x01;
                glm::vec3 row1(vec.x * vec.x, vec.x * vec.y, vec.x * vec.z);
                glm::vec3 row2(vec.x * vec.y, vec.y * vec.y, vec.y * vec.z);
                glm::vec3 row3(vec.x * vec.z, vec.y * vec.z, vec.z * vec.z);
                glm::mat3 hess(row1, row2, row3);
                float l = sqrt(glm::dot(x01, x01));
                hess = hess / l;
                glm::mat3 idm(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
                hess += (1 - spring.RestLength / l) * (idm - hess*l);
                hess *= -system.Stiffness;
                if (system.Fixed[p0] && system.Fixed[p1]) continue;
                forces[p0] += f;
                forces[p1] -= f;
                if (system.Fixed[p1]) {
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            coefficients.push_back(Eigen::Triplet<float>(3 * p0 + i, 3 * p0 + j, -hess[i][j]));
                        }
                    }
                } else if (system.Fixed[p0]) {
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            coefficients.push_back(Eigen::Triplet<float>(3 * p1 + i, 3 * p1 + j, -hess[i][j]));
                        }
                    }
                } else {
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            coefficients.push_back(Eigen::Triplet<float>(3 * p0 + i, 3 * p0 + j, -hess[i][j]));
                        }
                    }
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            coefficients.push_back(Eigen::Triplet<float>(3 * p1 + i, 3 * p1 + j, -hess[i][j]));
                        }
                    }
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            coefficients.push_back(Eigen::Triplet<float>(3 * p1 + i, 3 * p0 + j, hess[i][j]));
                        }
                    }
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            coefficients.push_back(Eigen::Triplet<float>(3 * p0 + i, 3 * p1 + j, hess[i][j]));
                        }
                    }
                }
                
            }

            Eigen::VectorXf b = Eigen::VectorXf::Zero(3 * n);
            for (int i = 0; i < n; i++) {
                for (int row = 0; row < 3; row++) {
                    coefficients.push_back(Eigen::Triplet(3 * i + row, 3 * i + row, system.Mass / ddt / ddt));
                }
                glm::vec3 y = system.Positions[i] + ddt * system.Velocities[i] + ddt * ddt * glm::vec3(0, -system.Gravity, 0);
                if (!system.Fixed[i]) {
                    b(3 * i) = ((system.Mass) / ddt / ddt * (system.Positions[i] - y) - forces[i]).x;
                    b(3 * i + 1) = ((system.Mass) / ddt / ddt * (system.Positions[i] - y) - forces[i]).y;
                    b(3 * i + 2) = ((system.Mass) / ddt / ddt * (system.Positions[i] - y) - forces[i]).z;
                }
            }
            M.setFromTriplets(coefficients.begin(), coefficients.end());
            auto solver = Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>(M);
            Eigen::VectorXf x = solver.solve(-b);
            for (int i = 0; i < n; i++) {
                system.Positions[i].x += x(3 * i);
                system.Positions[i].y += x(3 * i + 1);
                system.Positions[i].z += x(3 * i + 2);
                forces[i] = glm::vec3(0, -system.Gravity, 0);
            }
             for (auto const spring : system.Springs) {
                auto const p0  = spring.AdjIdx.first;
                auto const p1  = spring.AdjIdx.second;
                glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
                glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
                glm::vec3 const e01 = glm::normalize(x01);
                glm::vec3 f = (system.Stiffness * (glm::length(x01) - spring.RestLength) + system.Damping * glm::dot(v01, e01)) * e01;
                forces[p0] += f;
                forces[p1] -= f;
            }
            for (int i = 0; i < n; i++) {
                system.Velocities[i] += forces[i] * ddt / system.Mass;
            }
        }
    }
}
