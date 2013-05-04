#ifndef __PENDULUMS_H__
#define __PENDULUMS_H__

#include <vector>
#include <cmath>
#include <string>

namespace particles{
	const static double PI = 4*std::atan(1.0);
	const static double g = 0.0002;
	const static double r_threshold = 0.05;
	using std::sin;
	using std::cos;

	struct partState{
		double x, y;
		double vx, vy;
		partState()	: x(0.0), y(0.0), vx(0.0), vy(0.0){}
		partState(double _x, double _y, double _vx, double _vy)
			: x(_x), y(_y), vx(_vx), vy(_vy)
		{}
	};
	struct partSystem{
		std::vector<partState> pss;
		partSystem(int size){
			pss = std::vector<partState>(size);
		}
	};
	const partState operator+(const partState& lhs, const partState& rhs){
		partState ans;
		ans.x = lhs.x + rhs.x;
		ans.y = lhs.y + rhs.y;
		ans.vx = lhs.vx + rhs.vx;
		ans.vy = lhs.vy + rhs.vy;
		return ans;
	}
	const partState operator-(const partState& lhs, const partState& rhs){
		partState ans;
		ans.x = lhs.x - rhs.x;
		ans.y = lhs.y - rhs.y;
		ans.vx = lhs.vx - rhs.vx;
		ans.vy = lhs.vy - rhs.vy;
		return ans;
	}
	const partSystem operator+(const partSystem& lhs, const partSystem& rhs){
		partSystem ans(lhs.pss.size());
		for(unsigned int i = 0; i < lhs.pss.size(); i++){
			ans.pss.at(i) = lhs.pss.at(i) + rhs.pss.at(i);
		}
		return ans;
	}
	const partSystem operator-(const partSystem& lhs, const partSystem& rhs){
		partSystem ans(lhs.pss.size());
		for(unsigned int i = 0; i < lhs.pss.size(); i++){
			ans.pss.at(i) = lhs.pss.at(i) - rhs.pss.at(i);
		}
		return ans;
	}
	const partState operator*(const double lhs, const partState& rhs){
		partState ans;
		ans.x = lhs * rhs.x;
		ans.y = lhs * rhs.y;
		ans.vx = lhs * rhs.vx;
		ans.vy = lhs * rhs.vy;
		return ans;
	}
	const partSystem operator*(const double lhs, const partSystem& rhs){
		partSystem ans(rhs.pss.size());
		for(unsigned int i = 0; i < rhs.pss.size(); i++){
			ans.pss.at(i) = lhs * rhs.pss.at(i);
		}
		return ans;
	}
	class particleClass{
		int size;
		partSystem ps;
		std::vector<double> masses;
		int integrateMethod;
		double initEnergy;
	public:
		particleClass() : size(0), ps(0), integrateMethod(1), initEnergy(0.0){}
		particleClass(const std::vector<std::pair<double, partState> >& particles) : ps(particles.size()), integrateMethod(1){
			size = particles.size();
			double totalMass = 0.0;
			auto psit = ps.pss.begin();
			for(auto it = particles.begin(); it != particles.end(); ++it){
				totalMass += it->first;
				masses.push_back(it->first);
				*psit = it->second;
				++psit;
			}
			double vxDiff = mx() / totalMass;
			double vyDiff = my() / totalMass;
			for(auto it = ps.pss.begin(); it != ps.pss.end(); ++it){
				it->vx -= vxDiff;
				it->vy -= vyDiff;
			}
			initEnergy = energy();
		}
		~particleClass(){};
		void setIntegrate(int _integrateMethod){
			integrateMethod = _integrateMethod;
		}

		std::string getIntegrateName(){
			switch(integrateMethod){
			case 1:
				return "Euler Method";
			case 2:
				return "Heun's Method";
			case 3:
				return "Runge-Kutta Method";
			case 4:
				return "Velocity Velret Method";
			default:
				return "Unknown";
			}
		}
		std::vector<std::pair<double, double> > vertices(){
			std::vector<std::pair<double, double> > ans;
			for(auto it = ps.pss.begin(); it != ps.pss.end(); ++it){
				ans.push_back(std::pair<double, double>(it->x, it->y));
			}
			return ans;
		}
		void move(const double dt){
			switch(integrateMethod){
			case 1:
				move_euler(dt);
				break;
			case 2:
				move_heun(dt);
				break;
			case 3:
				move_runge_kutta(dt);
				break;
			case 4:
				move_velocity_velret(dt);
				break;
			}
		}
		void move_euler(const double dt){
			auto psDiff = calcAccel(ps);
			ps = ps + dt * psDiff;
		}
		void move_heun(const double dt){
			auto k1 = calcAccel(ps);
			auto k2 = calcAccel(ps + dt*k1);
			ps = ps + dt/2.0*(k1+k2);
		}
		void move_runge_kutta(const double dt){
			auto k1 = calcAccel(ps);
			auto k2 = calcAccel(ps + dt/2.0*k1);
			auto k3 = calcAccel(ps + dt/2.0*k2);
			auto k4 = calcAccel(ps + dt*k3);
			ps = ps + dt/6.0*(k1+2*k2+2*k3+k4);
		}
		void move_velocity_velret(const double dt){
			auto k1 = calcAccel(ps);
			auto p2 = ps;
			for(int i = 0; i < size; i++){
				p2.pss.at(i).x += dt*ps.pss.at(i).vx + dt*dt/2.0*k1.pss.at(i).vx;
				p2.pss.at(i).y += dt*ps.pss.at(i).vy + dt*dt/2.0*k1.pss.at(i).vy;
			}
			auto k2 = calcAccel(p2);
			for(int i = 0; i < size; i++){
				p2.pss.at(i).vx += dt/2.0*(k2.pss.at(i).vx + k1.pss.at(i).vx);
				p2.pss.at(i).vy += dt/2.0*(k2.pss.at(i).vy + k1.pss.at(i).vy);
			}
			ps = p2;
		}
		partSystem calcAccel(const partSystem& _ps){
			auto ans = calcForce(_ps);
			for(int i = 0; i < size; i++){
				ans.pss.at(i).x = _ps.pss.at(i).vx;
				ans.pss.at(i).y = _ps.pss.at(i).vy;
				ans.pss.at(i).vx /= masses.at(i);
				ans.pss.at(i).vy /= masses.at(i);
			}
			return ans;
		}
		partSystem calcForce(const partSystem& _ps){
			partSystem ans(size);
			for(int i = 0; i < size; i++){
				for(int j = i+1; j < size; j++){
					double dx = _ps.pss.at(i).x - _ps.pss.at(j).x, dy = _ps.pss.at(i).y - _ps.pss.at(j).y;
					double r2 = dx*dx + dy*dy;
					double r = std::sqrt(r2);
					if(r > r_threshold){
						//double r12_e3 = r2/0.0001;
						//r12_e3 = r12_e3*r12_e3*r12_e3*r12_e3*r12_e3*r12_e3;
						double f = -g*masses.at(i)*masses.at(j)/r;// + g2*masses.at(i)*masses.at(j)/r12_e3;
						double fx = f*dx/r, fy = f*dy/r;
						ans.pss.at(i).vx += fx;
						ans.pss.at(i).vy += fy;
						ans.pss.at(j).vx -= fx;
						ans.pss.at(j).vy -= fy;
					}
				}
				ans.pss.at(i).x = _ps.pss.at(i).vx;
				ans.pss.at(i).y = _ps.pss.at(i).vy;
			}
			return ans;
		}
		double energyDiff(){
			return energy()-initEnergy;
		}
		double energy(){
			return U()+T();
		}
		double U(){
			double ans = 0.0;
			for(int i = 0; i < size; i++){
				for(int j = i+1; j < size; j++){
					double dx = ps.pss.at(i).x - ps.pss.at(j).x, dy = ps.pss.at(i).y - ps.pss.at(j).y;
					double r2 = dx*dx + dy*dy;
					double r = std::sqrt(r2);
					r = std::max(r, r_threshold);
					ans -= g*masses.at(i)*masses.at(j)*(-std::log(r));
				}
			}
			return ans;
		}
		double T(){
			double ans = 0.0;
			for(int i = 0; i < size; i++){
				ans += 0.5*masses.at(i)*(ps.pss.at(i).vx*ps.pss.at(i).vx + ps.pss.at(i).vy*ps.pss.at(i).vy);
			}
			return ans;
		}
		double mx(){
			double ans = 0.0;
			for(int i = 0; i < size; i++){
				ans += masses.at(i)*ps.pss.at(i).vx;
			}
			return ans;
		}
		double my(){
			double ans = 0.0;
			for(int i = 0; i < size; i++){
				ans += masses.at(i)*ps.pss.at(i).vy;
			}
			return ans;
		}
	};
}

#endif
