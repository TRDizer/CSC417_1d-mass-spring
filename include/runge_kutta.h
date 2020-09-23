//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    Eigen::VectorXd f1, f2, f3, f4;
    Eigen::VectorXd q1 = q;
    Eigen::VectorXd qdot1 = qdot;
    force(f1, q1, qdot1);
    Eigen::VectorXd q2 = qdot1 * dt / 2 + q;
    Eigen::VectorXd qdot2 = f1 / mass * dt / 2 + qdot;
    force(f2, q2, qdot2);
    Eigen::VectorXd q3 = qdot2 * dt /2 + q;
    Eigen::VectorXd qdot3 = f2 / mass * dt / 2 + qdot;
    force(f3, q3, qdot3);
    Eigen::VectorXd q4 = qdot3 * dt + q;
    Eigen::VectorXd qdot4 = f3 / mass * dt + qdot;
    force(f4, q4, qdot4);
    q = q + dt * (qdot1 + 2*qdot2 + 2*qdot3 + qdot4) / 6;
    qdot = (f1 + 2*f2 + 2*f3 + f4) / 6 / mass * dt + qdot;
}