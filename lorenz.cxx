#include <iostream>
#include <cmath>

using namespace std;

void rk_4(double Y[], int N, double dt){

int a = 10;
int b = 28;
double c = 8.0/3;

double* k1 = new double [3];
double* k2 = new double [3];
double* k3 = new double [3];
double* k4 = new double [3];

Y[0] = 1;
Y[N] = 1;
Y[2*N] = 1;

for( int i=0 ; i<N-1 ; i++){
	// Calculation of k1
	k1[0] = a*(Y[N+i]-Y[i]);
	k1[1] = Y[i]*(b-Y[2*N+i])-Y[N+i];
	k1[2] = Y[i]*Y[N+i]-c*Y[2*N+i];
	// k1[0] = a*(y_n - x_n)
	// k1[1] = x_n*(b-z_n) - y_n
	// k1[2] = x_n*y_n - c*z_n
    //cout << "k1:" << "\t" << k1[0] << "\t" << k1[1] << "\t" << k1[2] << endl;

	// Calculation of k2
	k2[0] = a*( (Y[N+i]+0.5*dt*k1[0]) - (Y[i]+0.5*dt*k1[0]) );
	k2[1] = (Y[i]+0.5*dt*k1[0]) * (b - (Y[2*N+i]+0.5*dt*k1[0])) - (Y[N+i]+0.5*dt*k1[0]);
	k2[2] = (Y[i]+0.5*dt*k1[0]) * (Y[N+i]+0.5*dt*k1[0]) - c*(Y[2*N+i]+0.5*dt*k1[0]);
    //cout << "k2:" << "\t" << k2[0] << "\t" << k2[1] << "\t" << k2[2] << endl;

	// Calculation of k3
	k3[0] = a*( (Y[N+i]+0.5*dt*k2[0]) - (Y[i]+0.5*dt*k2[0]) );
	k3[1] = (Y[i]+0.5*dt*k2[0]) * (b - (Y[2*N+i]+0.5*dt*k2[0])) - (Y[N+i]+0.5*dt*k2[0]);
	k3[2] = (Y[i]+0.5*dt*k2[0]) * (Y[N+i]+0.5*dt*k2[0]) - c*(Y[2*N+i]+0.5*dt*k2[0]);
    //cout << "k3:" << "\t" << k3[0] << "\t" << k3[1] << "\t" << k3[2] << endl;

	// Calculation of k4
	k4[0] = a*( (Y[N+i]+dt*k3[0]) - (Y[i]+dt*k3[0]) );
	k4[1] = (Y[i]+dt*k3[0]) * (b - (Y[2*N+i]+dt*k3[0])) - (Y[N+i]+dt*k3[0]);
	k4[2] = (Y[i]+dt*k3[0]) * (Y[N+i]+dt*k3[0]) - c*(Y[2*N+i]+dt*k3[0]);
    //cout << "k4:" << "\t" << k4[0] << "\t" << k4[1] << "\t" << k4[2] << endl;

	// Calculation of y_(n+1) = y_(n) + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
	Y[i+1] = Y[i] + (1.0/6)*dt*( k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
	Y[N+i+1] = Y[N+i] + (1.0/6)*dt*( k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
	Y[2*N+i+1] = Y[2*N+i] + (1.0/6)*dt*( k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

	cout << i*dt << "\t" << Y[i] << "\t" << Y[N+i] << "\t" << Y[2*N+i] << endl;
    }
}

int main(){
double dt = 1.0/1600;
int T_end = 100;
int N = int(T_end/dt)+1;

double* Y = new double [3*N];

rk_4(Y, N, dt);

delete[] Y;

return 0;
}
